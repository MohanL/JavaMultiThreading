/*
    MST.java

    Euclidean minimum spanning tree implementation.

    Uses Dwyer's algorithm for Delaunay triangularization, then
    Kruskal's MST algorithm on the resulting mesh.  You need parallelize
    both stages.  For the triangulation stage, I recommend creating a
    new thread class similar to worker, to use in the divide-and-conquer
    step of triangulate().  For the tree state, I recommend letting worker
    threads find subtrees to merge in parallel, but forcing them to
    finalize the merges in order (after double-checking to make sure
    their subtrees haven't been changed by earlier finalizations).

    There are better ways to parallelize each of these steps, but
    they're quite a bit harder.

    Michael L. Scott, November 2015; based heavily on earlier
    incarnations of several programming projects, and on Delaunay mesh
    code written in 2007.
 */

/* Todo 
    2. triangulation divide and conquerthread done
    3. kruskal's
*/
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
import java.lang.*;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.List;

public class MST {
    private static int n = 50;              // default number of points
    private static long sd = 0;             // default random number seed
    private static int numThreads = 1;      // default 
    private static String fn = "";
    private static final int TIMING_ONLY    = 0;
    private static final int PRINT_EVENTS   = 1;
    private static final int SHOW_RESULT    = 2;
    private static final int FULL_ANIMATION = 3;
    private static int animate = TIMING_ONLY;       // default

    // Examine command-line arguments for alternative running modes.
    //
    private static void parseArgs(String[] args) throws Exception {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-a")) {
                if (++i >= args.length) {
                    System.err.print("Missing animation level\n");
                } else {
                    int an = -1;
                    try {
                        an = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (an >= TIMING_ONLY && an <= FULL_ANIMATION) {
                        animate = an;
                    } else {
                        System.err.printf("Invalid animation level: %s\n", args[i]);
                    }
                }
            } else if (args[i].equals("-n")) {
                if (++i >= args.length) {
                    System.err.print("Missing number of points\n");
                } else {
                    int np = -1;
                    try {
                        np = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (np > 0) {
                        n = np;
                    } else {
                        System.err.printf("Invalid number of points: %s\n", args[i]);
                    }
                }
            } else if (args[i].equals("-s")) {
                if (++i >= args.length) {
                    System.err.print("Missing seed\n");
                } else {
                    try {
                        sd = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) {
                        System.err.printf("Invalid seed: %s\n", args[i]);
                    }
                }
            } 
            else if (args[i].equals("-t")) {
                if (++i >= args.length) {
                    System.err.print("Missing number of threads\n");
                } else {
                    int nt = -1;
                    try {
                        nt = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (nt > 0) {
                        numThreads = nt;
                        // when I pass -t 1 as the cmd args, one thread got spawned
                    } else {
                        System.err.printf("Invalid number of threads: %s\n", args[i]);
                    }
                }
            }
            else if (args[i].equals("-c")) {
                if (++i >= args.length) {
                    System.err.print("Missing configuration files\n");
                } else {
                        fn = args[i];
                        FileReader fr = new FileReader(fn); 
                        BufferedReader br = new BufferedReader(fr); 
                        String s; 
                        while((s = br.readLine()) != null) 
                        { 
                             System.out.println("ckpt00");
                             System.out.println(s); 
                             String[] array = s.split(" ");
                             if ( array.length != 2 ) 
                             {
                                    System.out.printf("Invalid configuration file size = %s\n",array.length);
                             }
                             else
                             {
                                 System.out.printf("valid configuration file size = %s\n",array.length);
                             }
                             if (array[0].equals("-t"))
                             {
                                 int nt = Integer.parseInt(array[1]);
                                 if (nt > 0) {
                                    numThreads = nt;
                                 } else {
                                     System.err.printf("Invalid number of threads: %s\n", array[1]);
                                 }
                             }
                             else if (array[0].equals("-a"))
                             {
                                int an = -1;
                                try {
                                    an = Integer.parseInt(array[1]);
                                } catch (NumberFormatException e) { }
                                if (an >= TIMING_ONLY && an <= FULL_ANIMATION) {
                                    animate = an;
                                } else {
                                    System.err.printf("Invalid animation level: %s\n", array[1]);
                                }
                             }
                             else if (array[0].equals("-t")) 
                             {
                                int nt = -1;
                                try {
                                    nt = Integer.parseInt(args[i]);
                                } catch (NumberFormatException e) { }
                                if (nt > 0) {
                                    numThreads = nt;
                                    // when I pass -t 1 as the cmd args, one thread got spawned
                                } else {
                                    System.err.printf("Invalid number of threads: %s\n", array[1]);
                                } 
                             }                              
                             else if (array[0].equals("-s"))  
                             {
                                try {
                                    sd = Integer.parseInt(array[1]);
                                } catch (NumberFormatException e) {
                                    System.err.printf("Invalid seed: %s\n", array[1]);
                                }
                             }
                             else{
                                System.err.printf("doesn't match any, Invalid configuration file, element 0: %s  element 1: %s\n", array[0],array[1]);
                             }                                                     
                        } 
                        fr.close(); 
                }
            }
            else {
                System.err.printf("Unexpected argument: %s\n", args[i]);
            }
        }
    }

    // Initialize appropriate program components for specified animation mode.
    //
    // 
    private Surface build(RootPaneContainer pane, int an) {
        final Coordinator c = new Coordinator();
        Surface s = new Surface(n, sd, c, numThreads);
        Animation t = null;
        if (an == SHOW_RESULT || an == FULL_ANIMATION) {
            t = new Animation(s);
            new UI(c, s, t, sd, pane);
        }
        final Animation a = t;
        if (an == PRINT_EVENTS) {
            s.setHooks(
                new Surface.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("created   %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }},
                new Surface.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("destroyed %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }},
                new Surface.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("selected  %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }});
        } else if (an == FULL_ANIMATION) {
            Surface.EdgeRoutine er = new Surface.EdgeRoutine() {
                public void run(int x1, int y1, int x2, int y2, boolean dum)
                        throws Coordinator.KilledException {
                    c.hesitate();
                    a.repaint();        // graphics need to be re-rendered
                }};
            s.setHooks(er, er, er);
        }
        return s;
    }
    public static void main(String[] args) throws Exception {
        parseArgs(args);
        MST me = new MST();
        JFrame f = null;
        if (animate == SHOW_RESULT || animate == FULL_ANIMATION) {
            f = new JFrame("MST");
            f.addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    System.exit(0);
                }
            });
        } else {
            System.out.printf("%d points, seed %d\n", n, sd);
        }
        Surface s = me.build(f, animate);
        if (f != null) {
            f.pack();
            f.setVisible(true);
        } else {
            // Using terminal I/O rather than graphics.
            // Execute the guts of the run button handler method here.
            long startTime = new Date().getTime();
            long midTime = 0;
            try {
                s.DwyerSolve();
                midTime = new Date().getTime();
                //System.out.println("ckpt*");
                //System.out.println(s.edges.size());
                s.KruskalSolve();
            } catch(Coordinator.KilledException e) { }
            long endTime = new Date().getTime();
            System.out.printf("elapsed time: %.3f + %.3f = %.3f seconds\n",
                              (double) (midTime-startTime)/1000,
                              (double) (endTime-midTime)/1000,
                              (double) (endTime-startTime)/1000);
         }
      }
  }

// The Worker is the thread that does the actual work of finding a
// triangulation and MST (in the animated version -- main thread does it
// in the terminal I/O version).
//
   class Worker extends Thread {
    private final Surface s;
    private final Coordinator c;
    private final UI u;
    private final Animation a;

    // The run() method of a Java Thread is never invoked directly by
    // user code.  Rather, it is called by the Java runtime when user
    // code calls start().
    //
    // The run() method of a worker thread *must* begin by calling
    // c.register() and end by calling c.unregister().  These allow the
    // user interface (via the Coordinator) to pause and terminate
    // workers.  Note how the worker is set up to catch KilledException.
    // In the process of unwinding back to here we'll cleanly and
    // automatically release any monitor locks.  If you create new kinds
    // of workers (as part of a parallel solver), make sure they call
    // c.register() and c.unregister() properly.
    //
    public void run() {
        try {
            c.register();
            s.DwyerSolve();
            s.KruskalSolve();
            c.unregister();
        } catch(Coordinator.KilledException e) { }
        if (a != null) {
            // Tell the graphics event thread to unset the default
            // button when it gets a chance.  (Threads other than the
            // event thread cannot safely modify the GUI directly.)
            a.repaint();
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    u.getRootPane().setDefaultButton(null);
                    u.updateTime();
                }
            });
        }
    }

    // Constructor
    //
    public Worker(Surface S, Coordinator C, UI U, Animation A) {
        s = S;
        c = C;
        u = U;
        a = A;
      }
    }

// The Surface is the MST world, containing all the points and edges.
// contains triangulation and kruskal's

   class Surface{
    // Much of the logic in Dwyer's algorithm is parameterized by
    // directional (X or Y) and rotational (clockwise, counterclockwise)
    // orientation.  The following constants get plugged into the
    // parameter slots.

    private static final int xdim = 0;
    private static final int ydim = 1;
    private static final int ccw = 0;
    private static final int cw = 1;
    private static int numThreads = 1;
    private static int numAvaibleThreads = numThreads;
    private static int divcounter = 0;
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
        // newly created  ns.
    private final int n;  // number of points
    private final point points[];
        // main array of points, used for partitioning and rendering
    private final HashSet<point> pointHash;
        // Used to ensure that we never have two points directly on top of
        // each other.  See point.hashCode and point.equals below.
    public final SortedSet<edge> edges;
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
        // so this return value would change all the time
        // so I guess I can have a local copy in the thread class ? and compare later ? 

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
    //private synchronized void triangulate(int l, int r, int low0, int high0, int low1, int high1, int parity) throws Coordinator.KilledException {
    //     System.out.printf("%s %s %s %s %s %s %s\n",l, r, low0, high0, low1,high1,parity);

    private void triangulate(int l, int r, int low0, int high0, int low1, int high1, int parity) throws Coordinator.KilledException {
         System.out.printf("%s %s %s %s %s %s %s\n",l, r, low0, high0, low1,high1,parity);
         //System.out.println("ckpt1");
         //try{Thread.sleep(500);}
         //catch(InterruptedException e){System.out.println("ckpt5");}
         final int dim0;  final int dim1;
         final int dir0;  final int dir1;

         if (parity == 0) {
             dim0 = xdim;  dim1 = ydim;
             dir0 = ccw;   dir1 = cw;
         } else {
             dim0 = ydim;  dim1 = xdim;
             dir0 = cw;    dir1 = ccw;
         }
         //System.out.println("ckpt2");
         if (l == r) {
             System.out.println("case 1 return");
             return;
         }
         if (l == r-1) {
             new edge(points[l], points[r], null, null, dir1);
                 // direction doesn't matter in this case
             System.out.println("case 2 return");
             return;
         }
         if (l == r-2) {     // make single triangle
             edge e2 = new edge(points[l+1], points[r], null, null, dir1);
             edge e1 = new edge(points[l], points[l+1], null, e2, dir1);
             if (externAngle(points[l], points[l+1], points[r], dir0)) {
                 // new edge is dir0 of edge 1, dir1 of edge 2
                 new edge(points[l], points[r], e1, e2, dir0);
             } else {
                 // new edge is dir1 of edge 1, dir0 of edge 2
                 new edge(points[l], points[r], e1, e2, dir1);
             }
             System.out.println("case 3 return");
             return;
         }
         //System.out.println("ckpt3");
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
            System.out.println("Divide & conquer");
            
            if ((numThreads>= 2) && (numThreads == numAvaibleThreads)){
                System.out.println("multithreading round 1");
                divcounter++;
                t T1 = new t("1",1,l, i, low1, high1, low0, mid, 1-parity);
                t T2 = new t("2",2,j, r, low1, high1, mid, high0, 1-parity);
                T1.start();
                T2.start();
                try{
                    T2.join();
                    T1.join();
                }
                catch( InterruptedException e)
                {
                    System.out.println("Interrupted");
                }
                finally{
                System.out.println("ckpt5");
                System.out.printf("%s %s\n", numThreads,numAvaibleThreads); // check, available threads = num of total threadss
                }
            }
            else if (numAvaibleThreads >= 2)
            {    divcounter++;
                System.out.println("multithreading round unknowm" + numAvaibleThreads);
                t T1 = new t(String.valueOf(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,l, i, low1, high1, low0, mid, 1-parity);
                t T2 = new t(String.valueOf(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,j, r, low1, high1, mid, high0, 1-parity);
                T1.start();
                T2.start();
                try{
                    T1.join();
                    T2.join();
                }
                catch( InterruptedException e)
                {
                    System.out.println("Interrupted");
                }
                //triangulate(j, r, low1, high1, mid, high0, 1-parity);
            }
            else{
              divcounter++;
              System.out.println("single thread round " + divcounter);
              triangulate(l, i, low1, high1, low0, mid, 1-parity);
              triangulate(j, r, low1, high1, mid, high0, 1-parity);
            }
            
//**************************************************************************************************************************************************
            // volatile varibales
            // should I make this code atomic ?
            //synchronized(this){
            System.out.println("ckpt6");
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
            System.out.println("ckpt7");
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
            System.out.println("ckpt8");
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
            System.out.println("ckpt9");
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
            System.out.println("ckpt10");
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
                // another divide conquer ? 
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
            System.out.println("ckpt11");
//***************************************************************************************************************************************************
            //}
            }
    }

    // This is the actual MST calculation.
    // It relies on the fact that set "edges" is sorted by length, so
    // enumeration occurs shortest-to-longest.
    //

/*
The easiest strategy is probably to retain the serial iteration over edges, allow threads to identify the subtrees they want to merge (an O(log n) operation) in parallel, 
and then force the merges to actually complete in order 
(possibly starting over if they discover that a previous merge has changed which edges are in which subtrees). 
You may discover, however, that the condition synchronization for this strategy consumes more time than it saves; you’ll want to address this in your write-up. 
*/

// Not finished
// need to analyze how -t 2 is better than -t 10?
    public void KruskalSolve()
        throws Coordinator.KilledException 
        {
        System.out.println("Start Kruskals");
        int numTrees = n;
        //for (edge e : edges)
        r[] threadholder = new r[numThreads];
        r[] threadholder2 = new r[numThreads];
        point[] repold1 = new point[numThreads]; // old rep
        point[] repold2 = new point[numThreads]; // old rep
        point[] repnew1 = new point[numThreads]; // new rep
        point[] repnew2 = new point[numThreads]; // new rep
        point[] pointholder1 = new point[numThreads];
        point[] pointholder2 = new point[numThreads];
        int counter = 0;
        final Iterator<edge> it = edges.iterator();
        edge e = it.next();
        int flag2 = 1; // enable threading recycle
        int maxthreadsused = 0;
        do
        {   
            System.out.printf("%s %s\n",e.points[0].hashCode(),e.points[1].hashCode());
            point st1,st2;
            // this is the entry for multithreading
            if ((numThreads >= 2)&&(numAvaibleThreads == numThreads))
            {   
                // initilize 
                for(int i = 0;i+2<= numThreads;i+=2)
                {   System.out.printf("kruskal round %s\n",i);
                    // create multithreads here
                    threadholder[i] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,e,e.points[0]);
                    threadholder2[i] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,e,e.points[1]);
                    repold1[i] = threadholder[i].retrep();
                    repold2[i] = threadholder2[i].retrep();
                    if(i+2 == numThreads)
                    {
                        System.out.printf("maximum threads used %s\n",i+2);
                        maxthreadsused = i+2;
                        break;
                    }
                    else if (it.hasNext()) // assign more threads
                    {
                        e = it.next();
                    }
                }
                for(int i = 0;i <= maxthreadsused-2;i+=2)
                {
                    threadholder[i].start();
                    threadholder2[i].start();
                }   
                System.out.println("ckpt13");
                try
                {   
                    for(int i = 0;i <= maxthreadsused-2;i+=2)
                    {
                        threadholder[i].join();
                        threadholder2[i].join();
                    }   
                }
                catch(InterruptedException g)
                {
                    System.out.printf("error ckpt16\n");
                }
                System.out.println("ckpt14");
                for(int j = 0;j <= maxthreadsused-2;j+=2)
                {   
                    // conditional merge
                    repnew1[j] = threadholder[j].retrep();
                    repnew2[j] = threadholder2[j].retrep();

                    if((repold1[j] == repnew1[j]) && (repold2[j] == repnew2[j]))
                    {
                        pointholder1[j] = threadholder[j].ret();
                        pointholder2[j] = threadholder2[j].ret();

                        if (pointholder1[j] != pointholder2[j]) {
                            pointholder1[j].merge(pointholder2[j]);
                            threadholder[j].rete().addToMST();
                            System.out.printf("number of Trees in total %s\n", numTrees);
                            if (--numTrees == 1) {
                                break;
                            }
                        }
                    }
                   
                    else 
                    {
                        System.out.println("ckpt17");
                        edge g = threadholder[j].rete();
                        threadholder[j] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,g,g.points[0]);
                        threadholder2[j] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,g,g.points[1]);
                        threadholder[j].start();
                        threadholder2[j].start();
                        try
                        {   
                            threadholder[j].join();
                            threadholder2[j].join();
                        }   
                        catch(InterruptedException exc)
                        {
                            System.out.printf("error ckpt16\n");
                        }
                        pointholder1[j] = threadholder[j].ret();
                        pointholder2[j] = threadholder2[j].ret();

                        if (pointholder1[j] != pointholder2[j]) {
                            pointholder1[j].merge(pointholder2[j]);
                            threadholder[j].rete().addToMST();
                            System.out.printf("number of Trees in total %s\n", numTrees);
                            if (--numTrees == 1) {
                                break;
                            }
                        }

                    }
                    
                }
                System.out.println("ckpt18");
            }
            // in this case, some dead threads will be reused
            // currently there is no need to use this
            //***************************************************************************************************************************************************************************
            else if ((numAvaibleThreads >= 2) && (flag2 == 1))
            {   
                System.out.println("inside thread recycle");
                List<Integer> l = new ArrayList<Integer>();
                List<Integer> l2 = new ArrayList<Integer>();
                A:for(int i = 0;i+2<= maxthreadsused;i+=2)
                {   
                    boolean b1 = threadholder[i].isAlive();
                    boolean b2 = threadholder2[i].isAlive();
                    if((b1 == false) && (b2 == false))
                    {
                        threadholder[i] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,e,e.points[0]);
                        threadholder2[i] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,e,e.points[1]);
                        repold1[i] = threadholder[i].retrep();
                        repold2[i] = threadholder2[i].retrep();
                        l.add(i);
                        if (it.hasNext()) // assign more threads
                        {
                            e = it.next();
                        }
                        else
                        {
                            break A;
                        }
                    }
                    else
                    {
                        l2.add(i);
                    }
                }
                try
                {   
                    for(Integer i :l2)
                    {
                        threadholder[i].join();
                        threadholder2[i].join();
                    }   
                }
                catch(InterruptedException g)
                {
                    System.out.printf("error ckpt26\n");
                }

                for(Integer i:l)
                {
                    threadholder[i].start();
                    threadholder2[i].start();
                }   
                System.out.println("ckpt23");
                try
                {   
                    for(Integer i :l)
                    {
                        threadholder[i].join();
                        threadholder2[i].join();
                    }   
                }
                catch(InterruptedException g)
                {
                    System.out.printf("error ckpt26\n");
                }
                System.out.println("ckpt24");
                for(Integer j:l)
                {   
                    // conditional merge
                    repnew1[j] = threadholder[j].retrep();
                    repnew2[j] = threadholder2[j].retrep();

                    if((repold1[j] == repnew1[j]) && (repold2[j] == repnew2[j]))
                    {
                        System.out.println("ckpt25");
                        pointholder1[j] = threadholder[j].ret();
                        pointholder2[j] = threadholder2[j].ret();

                        if (pointholder1[j] != pointholder2[j]) {
                            pointholder1[j].merge(pointholder2[j]);
                            threadholder[j].rete().addToMST();
                            System.out.printf("number of Trees in total %s\n", numTrees);
                            if (--numTrees == 1) {
                                break;
                            }
                        }
                    }
                    else 
                    {
                        System.out.println("ckpt26");
                        edge g = threadholder[j].rete();
                        threadholder[j] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,g,g.points[0]);
                        threadholder2[j] = new r(Integer.toString(numThreads- numAvaibleThreads+1),numThreads-numAvaibleThreads+1,g,g.points[1]);
                        threadholder[j].start();
                        threadholder2[j].start();
                        try
                        {   
                            threadholder[j].join();
                            threadholder2[j].join();
                        }   
                        catch(InterruptedException exc)
                        {
                            System.out.printf("error ckpt16\n");
                        }
                        pointholder1[j] = threadholder[j].ret();
                        pointholder2[j] = threadholder2[j].ret();

                        if (pointholder1[j] != pointholder2[j]) {
                            pointholder1[j].merge(pointholder2[j]);
                            threadholder[j].rete().addToMST();
                            System.out.printf("number of Trees in total %s\n", numTrees);
                            if (--numTrees == 1) {
                                break;
                            }
                        }

                    }
                    System.out.println("ckpt27");
                }
                System.out.println("ckpt28");
            }

            //***************************************************************************************************************************************************************************
            else
            {   
                st1 = e.points[0].subtree();
                st2 = e.points[1].subtree();
                if (st1 != st2) 
                {
                    // This edge joins two previously separate subtrees.
                    st1.merge(st2);
                    e.addToMST();
                    System.out.printf("number of Trees in total %s\n", numTrees);
                    if (--numTrees == 1) {
                        break;
                    }
                }
            }
            if(it.hasNext())
                e = it.next();
            else
                break;
            System.out.println("ckpt29");
        }while(it.hasNext());
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
   // kat, You could try to implement a runnable class instead of thread class to understand the code better.
   class t extends Thread {
        private String threadName;
        private int threadNumber;
        private int locall;
        private int localr;
        private int locallow0;
        private int localhigh0;
        private int locallow1;
        private int localhigh1;
        private int localparity;
  
    t(String name,int nt,int tl, int tr, int tlow0, int thigh0, int tlow1, int thigh1, int tparity){
       numAvaibleThreads -= 1;
       threadName = name;
       threadNumber = nt;
       locall = tl;
       localr = tr;
       locallow0 = tlow0;
       localhigh0= thigh0;
       locallow1 = tlow1;
       localhigh1 = thigh1;
       localparity = tparity;
       System.out.println("Creating " +  threadName );
    }
    public void run() {
      System.out.println("Running " +  threadName );
      try {
            System.out.println("Inside thread " + threadName);
            Surface.this.triangulate(locall, localr, locallow0, localhigh0, locallow1, localhigh1, localparity);
      }
       catch (Coordinator.KilledException e){
           System.out.print("Thread" + threadName + " Coordinator.KilledException.");
           System.exit(1);
      }
     System.out.println("Thread " +  threadName + " exiting.");
     numAvaibleThreads+=1;
    }    
 }
//*******************************************************************************************************************
class r extends Thread {
    private String threadName;
    private int threadNumber;
    private point lp;
    private point rt;
    public edge le;
    public point rep;

    r(String name,int nt,edge e, point p){
       numAvaibleThreads -= 1;
       threadName = name;
       threadNumber = nt;
       lp = p;
       le = e;
       System.out.println("Creating " +  threadName);
    }

    public void run() {
       System.out.println("Running " +  threadName );
       System.out.println("Inside thread " + threadName);
       rt =  lp.subtree();
       System.out.println("Thread " +  threadName + " exiting.");
       numAvaibleThreads+=1;
    }
    // the point to make sure the representative point is consistent
    public point retrep()
    {
        rep = lp.subtree();
        return rep;
    }

    public point ret(){
        return rt;
    }

    public edge rete(){
        return le;
    }
  }
}
// Class Animation is the one really complicated sub-pane of the user interface.
// 
class Animation extends JPanel {
    private static final int width = 512;      // canvas dimensions
    private static final int height = 512;
    private static final int dotsize = 6;
    private static final int border = dotsize;
    private final Surface s;

    // The next two routines figure out where to render the dot
    // for a point, given the size of the animation panel and the spread
    // of x and y values among all points.
    //
    private int xPosition(int x) {
        return (int)
            (((double)x-(double)s.getMinx())*(double)width
                /((double)s.getMaxx()-(double)s.getMinx()))+border;
    }
    private int yPosition(int y) {
        return (int)
            (((double)s.getMaxy()-(double)y)*(double)height
                /((double)s.getMaxy()-(double)s.getMiny()))+border;
    }

    // The following method is called automatically by the graphics
    // system when it thinks the Animation canvas needs to be
    // re-displayed.  This can happen because code elsewhere in this
    // program called repaint(), or because of hiding/revealing or
    // open/close operations in the surrounding window system.
    //
    public void paintComponent(final Graphics g) {
        final Graphics2D g2 = (Graphics2D) g;

        super.paintComponent(g);    // clears panel
        s.forAllEdges(new Surface.EdgeRoutine() {
            public void run(int x1, int y1, int x2, int y2, boolean bold) {
                if (bold) {
                    g2.setPaint(Color.red);
                    g2.setStroke(new BasicStroke(4));
                } else {
                    g2.setPaint(Color.gray);
                    g2.setStroke(new BasicStroke(1));
                }
                g.drawLine(xPosition(x1), yPosition(y1),
                           xPosition(x2), yPosition(y2));
            }
        });
        s.forAllPoints(new Surface.PointRoutine() {
            public void run(int x, int y) {
                g2.setPaint(Color.blue);
                g.fillOval(xPosition(x)-dotsize/2, yPosition(y)-dotsize/2,
                           dotsize, dotsize);
            }
        });
    }

    // UI needs to call this routine when point locations have changed.
    //
    public void reset() {
        repaint();      // Tell graphics system to re-render.
    }

    // Constructor
    //
    public Animation(Surface S) {
        setPreferredSize(new Dimension(width+border*2, height+border*2));
        setBackground(Color.white);
        setForeground(Color.black);
        s = S;
        reset();
    }
}

// Class UI is the user interface.  It displays a Surface canvas above
// a row of buttons and a row of statistics.  Actions (event handlers)
// are defined for each of the buttons.  Depending on the state of the
// UI, either the "run" or the "pause" button is the default (highlighted in
// most window systems); it will often self-push if you hit carriage return.
//
class UI extends JPanel {
    private final Coordinator coordinator;
    private final Surface surface;
    private final Animation animation;

    private final JRootPane root;
    private static final int externalBorder = 6;

    private static final int stopped = 0;
    private static final int running = 1;
    private static final int paused = 2;

    private int state = stopped;
    private long elapsedTime = 0;
    private long startTime;

    private final JLabel time = new JLabel("time: 0");

    public void updateTime() {
        Date d = new Date();
        elapsedTime += (d.getTime() - startTime);
        time.setText(String.format("time: %d.%03d", elapsedTime/1000,
                                                    elapsedTime%1000));
    }

    // Constructor
    //
    public UI(Coordinator C, Surface S,
            Animation A, long SD, RootPaneContainer pane) {
        final UI ui = this;
        coordinator = C;
        surface = S;
        animation = A;

        final JPanel buttons = new JPanel();   // button panel
            final JButton runButton = new JButton("Run");
            final JButton pauseButton = new JButton("Pause");
            final JButton resetButton = new JButton("Reset");
            final JButton randomizeButton = new JButton("Randomize");
            final JButton quitButton = new JButton("Quit");

        final JPanel stats = new JPanel();   // statistics panel

        final JLabel seed = new JLabel("seed: " + SD + "   ");

        runButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (state == stopped) {
                    state = running;
                    root.setDefaultButton(pauseButton);
                    Worker w = new Worker(surface, coordinator,
                                          ui, animation);
                    Date d = new Date();
                    startTime = d.getTime();
                    w.start();
                } else if (state == paused) {
                    state = running;
                    root.setDefaultButton(pauseButton);
                    Date d = new Date();
                    startTime = d.getTime();
                    coordinator.toggle();
                }
            }
        });
        pauseButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (state == running) {
                    updateTime();
                    state = paused;
                    root.setDefaultButton(runButton);
                    coordinator.toggle();
                }
            }
        });
        resetButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                state = stopped;
                coordinator.stop();
                root.setDefaultButton(runButton);
                surface.reset();
                animation.reset();
                elapsedTime = 0;
                time.setText("time: 0");
            }
        });
        randomizeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                state = stopped;
                coordinator.stop();
                root.setDefaultButton(runButton);
                long v = surface.randomize();
                animation.reset();
                seed.setText("seed: " + v + "   ");
                elapsedTime = 0;
                time.setText("time: 0");
            }
        });
        quitButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.exit(0);
            }
        });

        // put the buttons into the button panel:
        buttons.setLayout(new FlowLayout());
        buttons.add(runButton);
        buttons.add(pauseButton);
        buttons.add(resetButton);
        buttons.add(randomizeButton);
        buttons.add(quitButton);

        // put the labels into the statistics panel:
        stats.add(seed);
        stats.add(time);

        // put the Surface canvas, the button panel, and the stats
        // label into the UI:
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setBorder(BorderFactory.createEmptyBorder(externalBorder,
            externalBorder, externalBorder, externalBorder));
        add(A);
        add(buttons);
        add(stats);

        // put the UI into the Frame
        pane.getContentPane().add(this);
        root = getRootPane();
        root.setDefaultButton(runButton);
    }
}
