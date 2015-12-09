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