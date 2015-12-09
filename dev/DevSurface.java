class Surface{
    // Much of the logic in Dwyer's algorithm is parameterized by
    // directional (X or Y) and rotational (clockwise, counterclockwise)
    // orientation.  The following constants get plugged into the
    // parameter slots.

    private static final int xdim = 0;
    private static final int ydim = 1;
    private static final int ccw = 0;
    private static final int cw = 1;
    private static final int numThreads = 1;
    private static int numAvaibleThreads = numThreads;
    private int minx;   // smallest x value among all points
    private int miny;   // smallest y value among all points
    private int maxx;   // largest x value among all points
    private int maxy;   // largest y value among all points
    public int getMinx() {return minx;}
    public int getMiny() {return miny;}
    public int getMaxx() {return maxx;}
    public int getMaxy() {return maxy;}

    // The following 7 fields are set by the Surface constructor.
    private final Coordinator coord;
        // Not needed at present, but will need to be passed to any
        // newly created workers.
    private final int n;  // number of points
    private final point points[];
        // main array of points, used for partitioning and rendering
    private final HashSet<point> pointHash;
        // Used to ensure that we never have two points directly on top of
        // each other.  See point.hashCode and point.equals below.
    private final SortedSet<edge> edges;
        // Used for rendering.  Ordering supports the KruskalSolve stage.
    private long sd = 0;
    private final Random prn;     // pseudo-random number generator

    // 3x3 determinant.  Called by 4x4 determinant.
    //
    private double det3(double a, double b, double c,
                        double d, double e, double f,
                        double g, double h, double i) {
        return a * (e*i - f*h)
             - b * (d*i - f*g)
             + c * (d*h - e*g);
    }

    // 4x4 determinant.  Called by encircled (below).
    //
    private double det4(double a, double b, double c, double d,
                        double e, double f, double g, double h,
                        double i, double j, double k, double l,
                        double m, double n, double o, double p) {
        return a * det3(f, g, h, j, k, l, n, o, p)
             - b * det3(e, g, h, i, k, l, m, o, p)
             + c * det3(e, f, h, i, j, l, m, n, p)
             - d * det3(e, f, g, i, j, k, m, n, o);
    }

    // swap points[i] and points[j]
    //
    private void swap(int i, int j) {
        point t = points[i];
        points[i] = points[j];
        points[j] = t;
    }

    // A point is a mesh/tree vertex.
    // It also serves, in the Kruskal stage, as a union-find set.
    // Its x and y coordinates are private and final.
    // Use the getCoord method to read their values.
    //
    private class point {
        private final int coordinates[] = new int[2];
        public edge firstEdge;

        // The following two fields are needed in the Kruskal stage (only).
        private point representative = null;    // equivalence set (subtree)
        private int subtreeSize = 1;

        // This is essentially the "find" of a union-find implementation.
        public point subtree() {
            point p = this;
            while (p.representative != null) p = p.representative;
            return p;
        }

        // And this is the "union".
        public void merge(point p) {
            // Make larger set the representative of the smaller set.
            assert representative == null;
            if (subtreeSize > p.subtreeSize) {
                subtreeSize += p.subtreeSize;
                p.representative = this;
            } else {
                p.subtreeSize += subtreeSize;
                representative = p;
            }
        }

        private int getCoord(int dim) {
            return coordinates[dim];
        }

        // Override Object.hashCode and Object.equals.
        // This way two points are equal (and hash to the same slot in
        // HashSet pointHash) if they have the same coordinates, even if they
        // are different objects.
        //
        public int hashCode() {
            return coordinates[xdim] ^ coordinates[ydim];
        }
        public boolean equals(Object o) {
            point p = (point) o;            // run-time type check
            return p.coordinates[xdim] == coordinates[xdim]
                && p.coordinates[ydim] == coordinates[ydim];
        }

        // Constructor
        //
        public point(int x, int y) {
            coordinates[xdim] = x;  coordinates[ydim] = y;
            // firstEdge == null
        }
    }

    // Signatures for things someone might want us to do with a point or
    // an edge (e.g., display it).
    // 
    public interface EdgeRoutine {
        public void run(int x1, int y1, int x2, int y2, boolean treeEdge)
            throws Coordinator.KilledException;
    }
    public interface PointRoutine{
        public void run(int x, int y);
    }

    public void forAllPoints(PointRoutine pr) {
        for (point p : points) {
            pr.run(p.getCoord(xdim), p.getCoord(ydim));
        }
    }
    public void forAllEdges(EdgeRoutine er) {
        for (edge e : edges) {
            try {
                er.run(e.points[0].getCoord(xdim),
                       e.points[0].getCoord(ydim),
                       e.points[1].getCoord(xdim),
                       e.points[1].getCoord(ydim), e.isMSTedge);
            } catch (Coordinator.KilledException f) { }
        }
    }

    // Routines to call when performing the specified operations:
    private static EdgeRoutine edgeCreateHook = null;
    private static EdgeRoutine edgeDestroyHook = null;
    private static EdgeRoutine edgeSelectHook = null;

    // The following is separate from the constructor to avoid a
    // circularity problem: when working in FULL_ANIMATION mode, the
    // Animation object needs a reference to the Surface object, and the
    // Surface object needs references to the hooks of the Animation object.
    //
    public void setHooks(EdgeRoutine ech, EdgeRoutine edh, EdgeRoutine esh) {
        edgeCreateHook = ech;
        edgeDestroyHook = edh;
        edgeSelectHook = esh;
    }

    // Edges encapsulate the bulk of the information about the triangulation.
    // Each edge contains references to its endpoints and to the next
    // edges clockwise and counterclockwise about those endpoints.
    //
    private class edge {
        public final point[] points = new point[2];
        public final edge[][] neighbors = new edge[2][2];
            // indexed first by edge end and then by rotational direction
        private boolean isMSTedge = false;
        public final double length;

        // Return index of point p within edge
        //
        public int indexOf(point p) {
            if (points[0] == p) return 0;
            if (points[1] == p) return 1;
            return -1;      // so I get an error if I use it
        }

        // utility routine for constructor
        //
        private void initializeEnd(point p, edge e, int end, int dir) {
            if (e == null) {
                neighbors[end][dir] = neighbors[end][1-dir] = this;
                p.firstEdge = this;
            } else {
                int i = e.indexOf(p);
                neighbors[end][1-dir] = e;
                neighbors[end][dir] = e.neighbors[i][dir];
                e.neighbors[i][dir] = this;
                i = neighbors[end][dir].indexOf(p);
                neighbors[end][dir].neighbors[i][1-dir] = this;
            }
        }

        // Constructor: connect points A and B, inserting dir (CW or CCW)
        // of edge Ea at the A end and 1-dir of edge Eb at the B end.
        // Either or both of Ea and Eb may be null.
        //
        public edge(point A, point B, edge Ea, edge Eb, int dir)
                throws Coordinator.KilledException {
            points[0] = A;  points[1] = B;
            double dx = (double) A.getCoord(xdim) - (double) B.getCoord(xdim);
            double dy = (double) A.getCoord(ydim) - (double) B.getCoord(ydim);
            length = Math.sqrt(dx * dx + dy * dy);

            initializeEnd(A, Ea, 0, dir);
            initializeEnd(B, Eb, 1, 1-dir);

            edges.add(this);
            if (edgeCreateHook != null)
                edgeCreateHook.run(points[0].getCoord(xdim),
                                   points[0].getCoord(ydim),
                                   points[1].getCoord(xdim),
                                   points[1].getCoord(ydim), false);
        }

        // Destructor: take self out of edges, point edge lists.
        // Should only be called when flipping an edge, so destroyed
        // edge should have neighbors at both ends.
        //
        public void destroy() throws Coordinator.KilledException {
            edges.remove(this);
            for (int i = 0; i < 2; i++) {
                int cw_index = neighbors[i][cw].indexOf(points[i]);
                int ccw_index = neighbors[i][ccw].indexOf(points[i]);
                neighbors[i][cw].neighbors[cw_index][ccw] = neighbors[i][ccw];
                neighbors[i][ccw].neighbors[ccw_index][cw] = neighbors[i][cw];
                if (points[i].firstEdge == this)
                    points[i].firstEdge = neighbors[i][ccw];
            }
            if (edgeDestroyHook != null)
                edgeDestroyHook.run(points[0].getCoord(xdim),
                                    points[0].getCoord(ydim),
                                    points[1].getCoord(xdim),
                                    points[1].getCoord(ydim), false);
        }

        // Assume edges are unique.
        // Override Object.equals to make it consistent with
        // edgeComp.compare below.
        //
        public boolean equals(Object o) {
            return this == o;
        }

        // Label this edge as an MST edge.
        //
        public void addToMST() throws Coordinator.KilledException {
            isMSTedge = true;
            if (edgeSelectHook != null)
                edgeSelectHook.run(points[0].getCoord(xdim),
                                   points[0].getCoord(ydim),
                                   points[1].getCoord(xdim),
                                   points[1].getCoord(ydim), false);
        }
    }

    // To support ordered set of edges.  Return 0 _only_ if two
    // arguments are the _same_ edge (this is necessary for unique
    // membership in set).  Otherwise order based on length.  If lengths
    // are the same, order by coordinates.
    //
                                                            
    public static class edgeComp implements Comparator<edge> {

        public int compare(edge e1, edge e2) {
            if (e1.equals(e2)) return 0;
            if (e1.length < e2.length) return -1;
            if (e1.length > e2.length) return 1;
            int e1xmin = e1.points[0].getCoord(xdim)
                            < e1.points[1].getCoord(xdim) ?
                                e1.points[0].getCoord(xdim) :
                                e1.points[1].getCoord(xdim);
            int e2xmin = e2.points[0].getCoord(xdim)
                            < e2.points[1].getCoord(xdim) ?
                                e2.points[0].getCoord(xdim) :
                                e2.points[1].getCoord(xdim);
            if (e1xmin < e2xmin) return -1;
            if (e1xmin > e2xmin) return 1;
            int e1ymin = e1.points[0].getCoord(ydim)
                            < e1.points[1].getCoord(ydim) ?
                                e1.points[0].getCoord(ydim) :
                                e1.points[1].getCoord(ydim);
            int e2ymin = e2.points[0].getCoord(ydim)
                            < e2.points[1].getCoord(ydim) ?
                                e2.points[0].getCoord(ydim) :
                                e2.points[1].getCoord(ydim);
            if (e1ymin < e2ymin) return -1;
            // if (e1ymin > e2ymin)
                return 1;
            // no other options; endpoints have to be distinct
        }
    }

    // Called by the UI when it wants to reset with a new seed.
    //
     public long randomize() {
         sd++;
         reset();
         return sd;
     }

     // Called by the UI when it wants to start over.
     //
     public void reset() {
         prn.setSeed(sd);
         minx = Integer.MAX_VALUE;
         miny = Integer.MAX_VALUE;
         maxx = Integer.MIN_VALUE;
         maxy = Integer.MIN_VALUE;
         pointHash.clear();      // empty out the set of points
         for (int i = 0; i < n; i++) {
             point p;
             int x;
             int y;
             do {
                 x = prn.nextInt();
                 y = prn.nextInt();
                 p = new point(x, y);
             } while (pointHash.contains(p));
             pointHash.add(p);
             if (x < minx) minx = x;
             if (y < miny) miny = y;
             if (x > maxx) maxx = x;
             if (y > maxy) maxy = y;
             points[i] = p;
         }
         edges.clear();      // empty out the set of edges
     }

     // If A, B, and C are on a circle, in counter-clockwise order, then
     // D lies within that circle iff the following determinant is positive:
     //
     // | Ax  Ay  Ax^2+Ay^2  1 |
     // | Bx  By  Bx^2+By^2  1 |
     // | Cx  Cy  Cx^2+Cy^2  1 |
     // | Dx  Dy  Dx^2+Dy^2  1 |
     //
     private boolean encircled(point A, point B, point C, point D, int dir) {
         if (dir == cw) {
             point t = A;  A = C;  C = t;
         }
         double Ax = A.getCoord(xdim);   double Ay = A.getCoord(ydim);
         double Bx = B.getCoord(xdim);   double By = B.getCoord(ydim);
         double Cx = C.getCoord(xdim);   double Cy = C.getCoord(ydim);
         double Dx = D.getCoord(xdim);   double Dy = D.getCoord(ydim);

         return det4(Ax, Ay, (Ax*Ax + Ay*Ay), 1,
                     Bx, By, (Bx*Bx + By*By), 1,
                     Cx, Cy, (Cx*Cx + Cy*Cy), 1,
                     Dx, Dy, (Dx*Dx + Dy*Dy), 1) > 0;
     }

     // Is angle from p1 to p2 to p3, in direction dir
     // around p2, greater than or equal to 180 degrees?
     //
     private boolean externAngle(point p1, point p2, point p3, int dir) {
         if (dir == cw) {
             point t = p1;  p1 = p3;  p3 = t;
         }
         int x1 = p1.getCoord(xdim);     int y1 = p1.getCoord(ydim);
         int x2 = p2.getCoord(xdim);     int y2 = p2.getCoord(ydim);
         int x3 = p3.getCoord(xdim);     int y3 = p3.getCoord(ydim);

         if (x1 == x2) {                     // first segment vertical
             if (y1 > y2) {                  // points down
                 return (x3 >= x2);
             } else {                        // points up
                 return (x3 <= x2);
             }
         } else {
             double m = (((double) y2) - y1) / (((double) x2) - x1);
                 // slope of first segment
             if (x1 > x2) {      // points left
                 return (y3 <= m * (((double) x3) - x1) + y1);
                 // p3 below line
             } else {            // points right
                 return (y3 >= m * (((double) x3) - x1) + y1);
                 // p3 above line
             }
         }
     }

     // Divide points[l..r] into two partitions.  Solve recursively, then
     // stitch back together.  Dim0 values range from [low0..high0].
     // Dim1 values range from [low1..high1].  We partition based on dim0.
     // Base case when 1, 2, or 3 points.
     //
     // As suggested by Dwyer, we swap axes and rotational directions
     // at successive levels of recursion, to minimize the number of long
     // edges that are likely to be broken when stitching.
     //
     
    private void triangulate(int l, int r, int low0, int high0, int low1, int high1, int parity) throws Coordinator.KilledException {

         final int dim0;  final int dim1;
         final int dir0;  final int dir1;

         if (parity == 0) {
             dim0 = xdim;  dim1 = ydim;
             dir0 = ccw;   dir1 = cw;
         } else {
             dim0 = ydim;  dim1 = xdim;
             dir0 = cw;    dir1 = ccw;
         }

         if (l == r) {
             return;
         }
         if (l == r-1) {
             new edge(points[l], points[r], null, null, dir1);
                 // direction doesn't matter in this case
             return;
         }
         if (l == r-2) {     // make single triangle
             // how do I support more 2 threads here ? 
             // Solve recursively, do I need to extends edge class ?
             edge e2 = new edge(points[l+1], points[r], null, null, dir1);
             edge e1 = new edge(points[l], points[l+1], null, e2, dir1);
             if (externAngle(points[l], points[l+1], points[r], dir0)) {
                 // new edge is dir0 of edge 1, dir1 of edge 2
                 new edge(points[l], points[r], e1, e2, dir0);
             } else {
                 // new edge is dir1 of edge 1, dir0 of edge 2
                 new edge(points[l], points[r], e1, e2, dir1);
             }
             return;
         }

         // At this point we know we're not a base case; have to subdivide.

         int mid = low0/2 + high0/2;
         int i = l;  int j = r;

         point lp = points[l];          // rightmost point in left half;
         int lp0 = Integer.MIN_VALUE;   // X coord of lp
         point rp = points[r];          // leftmost point in right half;
         int rp0 = Integer.MAX_VALUE;   // X coord of rp

         while (true) {
             // invariants: [i..j] are unexamined;
             // [l..i) are all <= mid; (j..r] are all > mid.

             int i0 = 0;  int j0 = 0;

             while (i < j) {
                 i0 = points[i].getCoord(dim0);
                 if (i0 > mid) {     // belongs in right half
                     if (i0 < rp0) {
                         rp0 = i0;  rp = points[i];
                     }
                     break;
                 } else {
                     if (i0 > lp0) {
                         lp0 = i0;  lp = points[i];
                     }
                 }
                 i++;
             }

             while (i < j) {
                 j0 = points[j].getCoord(dim0);
                 if (j0 <= mid) {    // belongs in left half
                     if (j0 > lp0) {
                         lp0 = j0;  lp = points[j];
                     }
                     break;
                 } else {
                     if (j0 < rp0) {
                         rp0 = j0;  rp = points[j];
                     }
                 }
                 j--;
             }

             // at this point either i == j == only unexamined element
             // or i < j (found elements that need to be swapped)
             // or i = j+1 (and all elements are in order)
             if (i == j) {
                 i0 = points[i].getCoord(dim0);
                 if (i0 > mid) {
                     // give border element to right half
                     if (i0 < rp0) {
                         rp0 = i0;  rp = points[i];
                     }
                     i--;
                 } else {
                     // give border element to left half
                     if (i0 > lp0) {
                         lp0 = i0;  lp = points[i];
                     }
                     j++;
                 }
                 break;
             }
             if (i > j) {
                 i--;  j++;  break;
             }
             swap(i, j);
             i++;  j--;
         }
         // Now [l..i] is the left partition and [j..r] is the right.
         // Either may be empty.

         if (i < l) {
             // empty left half
             triangulate(j, r, low1, high1, mid, high0, 1-parity);
         } else if (j > r) {
             // empty right half
             triangulate(l, i, low1, high1, low0, mid, 1-parity);
         } else {
            // divide and conquer
            if ((numThreads> 2) && (numThreads == numAvaibleThreads)){
                t T1 = new t("1",1,l, i, low1, high1, low0, mid, 1-parity);
                T1.start();
                numAvaibleThreads -= 1;
                t T2 = new t("2",2,l, i, low1, high1, low0, mid, 1-parity);
                T2.start();
                numAvaibleThreads -= 1;
            }
            else if (numAvaibleThreads >= 2)
            {   
                numAvaibleThreads -= 1;
                t T1 = new t(String.valueOf(numThreads- numAvaibleThreads),numThreads-numAvaibleThreads,l, i, low1, high1, low0, mid, 1-parity);
                T1.start();
                numAvaibleThreads -= 1;
                t T2 = new t(String.valueOf(numThreads- numAvaibleThreads),numThreads-numAvaibleThreads,l, i, low1, high1, low0, mid, 1-parity);
                T2.start();
            }
            else{
              triangulate(l, i, low1, high1, low0, mid, 1-parity);
              triangulate(j, r, low1, high1, mid, high0, 1-parity);
            }
             // prepare to stitch meshes together up the middle:
             class side {
                 public point p;     // working point
                 public edge a;      // above p
                 public edge b;      // below p
                 public point ap;    // at far end of a
                 public point bp;    // at far end of b
                 public int ai;      // index of p within a
                 public int bi;      // index of p within b
             };
             side left = new side();
             side right = new side();
             left.p = lp;
             right.p = rp;

            // Rotate around extreme point to find edges adjacent to Y
            // axis.  This class is basically a hack to get around the
            // lack of nested subroutines in Java.  We invoke its run
            // method twice below.
            //
            class rotateClass {
                void run(side s, int dir) {
                    // rotate around s.p to find edges adjacent to Y axis
                    if (s.p.firstEdge != null) {
                        s.a = s.p.firstEdge;
                        s.ai = s.a.indexOf(s.p);
                        s.ap = s.a.points[1-s.ai];
                        if (s.a.neighbors[s.ai][dir] == s.a) {
                            // only one incident edge on the right
                            s.b = s.a;
                            s.bi = s.ai;
                            s.bp = s.ap;
                        } else {
                            // >= 2 incident edges on the right;
                            // need to find correct ones
                            while (true) {
                                s.b = s.a.neighbors[s.ai][dir];
                                s.bi = s.b.indexOf(s.p);
                                s.bp = s.b.points[1-s.bi];
                                if (externAngle(s.ap, s.p, s.bp, dir)) break;
                                s.a = s.b;
                                s.ai = s.bi;
                                s.ap = s.bp;
                            }
                        }
                    }
                }
            }
            rotateClass rotate = new rotateClass();
            rotate.run(left, dir1);
            rotate.run(right, dir0);

            // Find endpoint of bottom edge of seam, by moving around border
            // as far as possible without going around a corner.  This, too,
            // is basically a nested subroutine.
            //
            class findBottomClass {
                boolean move(side s, int dir, point o) {
                    boolean progress = false;
                    if (s.b != null) {
                        while (!externAngle(s.bp, s.p, o, 1-dir)) {
                            // move s.p in direction dir
                            progress = true;
                            s.a = s.b;
                            s.ai = 1-s.bi;
                            s.ap = s.p;
                            s.p = s.b.points[1-s.bi];
                            s.b = s.b.neighbors[1-s.bi][dir];
                            s.bi = s.b.indexOf(s.p);
                            s.bp = s.b.points[1-s.bi];
                        }
                    }
                    return progress;
                }
            }
            findBottomClass findBottom = new findBottomClass();
            do {} while (findBottom.move(left, dir1, right.p)
                      || findBottom.move(right, dir0, left.p));

            // create bottom edge:
            edge base = new edge(left.p, right.p,
                                 left.a == null ? left.b : left.a,
                                 right.a == null ? right.b : right.a,
                                 dir1);
            final edge bottom = base;
            if (left.a == null) left.a = bottom;
                // left region is a singleton
            if (right.a == null) right.a = bottom;
                // right region is a singleton

            // Work up the seam creating new edges and deleting old
            // edges where necessary.  Note that {left,right}.{b,bi,bp}
            // are no longer needed.

            while (true) {

                // Find candidate endpoint.  Yet another nested subroutine.
                //
                class findCandidateClass {
                    point call(side s, int dir, edge base, point o)
                            throws Coordinator.KilledException {
                            // o is at far end of base
                        if (s.a == bottom) {
                            // region is a singleton
                            return null;
                        }
                        point c = s.a.points[1-s.ai];
                        if (externAngle(o, s.p, c, dir)) {
                            // no more candidates
                            return null;
                        }
                        while (true) {
                            edge na = s.a.neighbors[s.ai][dir];
                                // next edge into region
                            if (na == base) {
                                // wrapped all the way around
                                return c;
                            }
                            int nai = na.indexOf(s.p);
                            point nc = na.points[1-nai];
                                // next potential candidate
                            if (encircled(o, c, s.p, nc, dir)) {
                                // have to break an edge
                                s.a.destroy();
                                s.a = na;
                                s.ai = nai;
                                c = nc;
                            } else return c;
                        }
                    }
                }
                findCandidateClass findCandidate = new findCandidateClass();
                point lc = findCandidate.call(left, dir0, bottom, right.p);
                point rc = findCandidate.call(right, dir1, bottom, left.p);

                if (lc == null && rc == null) {
                    // no more candidates
                    break;
                }
                // Choose between candidates:
                if (lc != null && rc != null &&
                        encircled (right.p, lc, left.p, rc, dir0)) {
                    // Left candidate won't work; circumcircle contains
                    // right candidate.
                    lc = null;
                }
                // Now we know one candidate is null and the other is not.
                if (lc == null) {
                    // use right candidate
                    right.a = right.a.neighbors[1-right.ai][dir1];
                    right.ai = right.a.indexOf(rc);
                    right.ap = right.a.points[1-right.ai];
                    right.p = rc;
                    base = new edge(left.p, rc, left.a, right.a, dir1);
                } else {
                    // use left candidate
                    left.a = left.a.neighbors[1-left.ai][dir0];
                    left.ai = left.a.indexOf(lc);
                    left.ap = left.a.points[1-left.ai];
                    left.p = lc;
                    base = new edge(lc, right.p, left.a, right.a, dir1);
                }
            }
        }
    }

    // This is the actual MST calculation.
    // It relies on the fact that set "edges" is sorted by length, so
    // enumeration occurs shortest-to-longest.
    //
    public void KruskalSolve()
        throws Coordinator.KilledException {
        int numTrees = n;
        for (edge e : edges) {
            point st1 = e.points[0].subtree();
            point st2 = e.points[1].subtree();
            if (st1 != st2) {
                // This edge joins two previously separate subtrees.
                st1.merge(st2);
                e.addToMST();
                if (--numTrees == 1) break;
            }
        }
    }

    // This is a wrapper for the root call to triangulate().
    //
    public void DwyerSolve() throws Coordinator.KilledException {
        triangulate(0, n-1, minx, maxx, miny, maxy, 0);
    }

    // Constructor
    //
    public Surface(int N, long SD, Coordinator C, int nt) {
        n = N;
        sd = SD;
        coord = C;
        numThreads = nt;
        numAvaibleThreads = nt;
        points = new point[n];
        edges = new ConcurrentSkipListSet<edge>(new edgeComp());
            // Supports safe concurrent access by worker and graphics threads,
            // and as a SortedSet it keeps the edges in order by length.
        pointHash = new HashSet<point>(n);

        prn = new Random();
        reset();
    }
//*******************************************************************************************************************
   class t extends Thread {
        private Thread t;
        private String threadName;
        private int threadNumber;
        private int l;
        private int r;
        private int low0;
        private int high0;
        private int low1;
        private int high1;
        private int parity;
    //int l, int r, int low0, int high0, int low1, int high1, int parity) throws Coordinator.KilledException 
    t(String name,int nt,int l, int r, int low0, int high0, int low1, int high1, int parity){
       threadName = name;
       threadNumber = nt;
       l = this.l;
       r = this.r;
       low0 = this.low0;
       high0= this.high0;
       low1 = this.low1;
       high1 = this.high1;
       parity = this.parity;
       System.out.println("Creating " +  threadName );
    }

    public void run() {
      System.out.println("Running " +  threadName );
      try {
            triangulate(l, i, low1, high1, low0, mid, 1-parity);
            //Thread.sleep(50);
     } catch (InterruptedException e) {
         System.out.println("Thread " +  threadName + " interrupted.");
     }
     System.out.println("Thread " +  threadName + " exiting.");
    }
   
    public void start (){
      System.out.println("Starting " +  threadName );
      if (t == null)
      {
         t = new Thread (this, threadName);
         t.start ();
      }
     }
    
 }
//*******************************************************************************************************************
    class R implements Runnable {
        private Thread t;
        private String threadName;
        private int threadNumber;

    t( String name,int nt){
       threadName = name;
       threadNumber = nt
       System.out.println("Creating " +  threadName );
    }

    public void run() {
      System.out.println("Running " +  threadName );
      try {
            // call triangulate
            Thread.sleep(50);
     } catch (InterruptedException e) {
         System.out.println("Thread " +  threadName + " interrupted.");
     }
     System.out.println("Thread " +  threadName + " exiting.");
    }
   
    public void start ()
    {
      System.out.println("Starting " +  threadName );
      if (t == null)
      {
         t = new Thread (this, threadName);
         t.start ();
      }
    }
  }
}