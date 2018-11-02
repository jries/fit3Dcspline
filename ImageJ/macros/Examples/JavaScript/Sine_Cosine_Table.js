// This script displays a sine/cosine table in a window.

   rt = new ResultsTable();
   for (n=0; n<=2*Math.PI; n += 0.1) {
      rt.incrementCounter();
      rt.addValue("n", n);
      rt.addValue("Sine(n)", Math.sin(n));
      rt.addValue("Cos(n)", Math.cos(n));
   }
   rt.showRowNumbers(false);
   rt.show("Sine/Cosine Table");
