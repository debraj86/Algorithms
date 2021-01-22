//block searched stochastic simulation algorithm (BlSSSA) developed by Debraj Ghosh
import java.lang.*;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.io.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Scanner;
import java.util.Random;
class BlSSSA
{
 Scanner sc = new Scanner(System.in);
  int M = sc.nextInt(), N = sc.nextInt();
     
   int trials=0; 
   int[][] a = new int [M][N]; int up=0; 
   double[] k = new double [N], /*AND = new double [M],*/ sigma = new double [M];  int[] conc = new int[M], conc_L = new int[M], conc_H = new int[M]; double fluc;
   int[] no_spec_per_reac = new int[N], no_input_spec_per_reac = new int[N], /*no_depend_per_reac = new int[N],*/ interacttran = new int [M], alreadyused = new int[N]/*, incrment_U_3 = new int[M]*/;
    int [][] spec_per_reac = new int [N][], input_spec_per_reac = new int [N][], change_input_conc = new int[N][], change_conc = new int[N][], /*depend_per_reac = new int[N][],*/ lookuptable = new int[M][]; 
    double[][] pie = new double[M][], pie_H = new double[M][], pie_L = new double[M][];  /*int[][][] U_3 = new int[M][][]; */
    double[] prop = new double[N];  double totalprop = 0, deltotalprop=0;  int ival=0, jval=0; double transitionfiretime, u, totaltransitionfiringtime=0;
  
Random rand = new Random();
public void slow_update(int i)
{
  
  if (conc[i] <= 0)
  {
    conc_L[i]=0;
    conc_H[i]=0;
  }
  else if(conc[i] < 25)
  {
    conc_H[i]=conc[i]+4;
    if((conc[i]-4) > 0)
      conc_L[i]=conc[i]-4;
    else
      conc_L[i]=0;
  }
  
  
  else
  {
    //conc_L[i]= conc[i]*9/10;
    //conc_H[i]= conc[i]*11/10;
    conc_L[i]= (int)((double)conc[i]*(1.0-fluc));
    conc_H[i]= (int)((double)conc[i]*(1.0+fluc));
  }
  
}

public  int firingtran(double[] intervalprop, int[] interval, int intervals)
{
   
  double sam, selector=0;  int xyz = 0; 
  
  int accepted = 0; u=1;
  while(accepted != 1)
{ trials = trials + 1; 
  
  sam = totalprop*Math.random();
  selector=0; 
  
  for(int i=0; i<intervals; i++)
  {
    selector = selector + intervalprop[i]; 
    if(selector > sam )
    {
    selector = selector - intervalprop[i]; 
    xyz = interval[i]; 
    break; 
    }
  }
  for (int i=xyz;i<M;i++)
  {
    selector = selector + sigma[i];
    if(selector > sam)
    {
      
      selector = selector - sigma[i];
      ival = i; 
      break; 
    }
  }
  
  for(int j = 0; j < interacttran[ival]; j++)
  {
    selector = selector + pie_H[ival][j];
    
    if (selector > sam)
    {
      jval = j; 
     break; 
    }
  }
  double r1 = rand.nextDouble();
  if((pie_L[ival][jval] != 0) && (r1 <= pie_L[ival][jval]/pie_H[ival][jval]))
  {
  accepted = 1; 
  break; 
  }
  else 
  {
  int x = lookuptable[ival][jval];
  pie[ival][jval] = conc[ival]*k[x]; 
  for(int k = 0; k < no_input_spec_per_reac[x]; k++)
  {
	 if ( input_spec_per_reac[x][k] != ival )
	 {
	   pie[ival][jval] = pie[ival][jval]*conc[input_spec_per_reac[x][k]];
	   
	 }
	 if( input_spec_per_reac[x][k] == ival && change_input_conc[x][k] == -2)
	 {
	   pie[ival][jval] = pie[ival][jval]*(conc[input_spec_per_reac[x][k]]-1)*0.5;
	   
	 }
  }
  if((pie[ival][jval] != 0) && (r1 <= pie[ival][jval]/pie_H[ival][jval]))
  {
  accepted = 1; 
  break; 
  }
  }
  u=u*Math.random();
  
}

  return lookuptable[ival][jval];
}
  
  public  void exec(int selrxn)
  {
    for(int i=0; i< no_spec_per_reac[selrxn]; i++)
    {
      int x = spec_per_reac[selrxn][i]; 
      conc[x]+= change_conc[selrxn][i]; 
      if (conc[x] < 0)
	conc[x] = 0; 
      
    }
   
    
  }

  public  void calcpropensity()
  {
    totalprop = 0; int j1 = 0;
    for(int i=0;i<N;i++)
      prop[i]=1;
    for(int i=0; i<M; i++)
      interacttran[i]=0; 
    for(int i=0; i<N; i++)
      alreadyused[i]=0; 
    for(int i = 0; i< M; i++)
    {
      for(int j=0; j<N; j++)
      {  
	
	if(a[i][j] < 0 && alreadyused[j] == 0) 
	{
	  interacttran[i]++; 	//The lookup table L for storing the groups of reactions indexed by the species 
	  
	  alreadyused[j] = 1;
	}
      }
    }
    for(int i=0; i<N; i++)
      alreadyused[i]=0;
  /*    System.out.println(); 
    for(int i=0; i<M; i++)
      System.out.print(" "+interacttran[i]);*/
    for(int i = 0; i< M; i++)
    {
      lookuptable[i] = new int[interacttran[i]];  
    }  
    for(int i = 0; i< M; i++)
    {
      pie[i] = new double[interacttran[i]];
      pie_H[i] = new double[interacttran[i]];         // lambda vector
      pie_L[i] = new double[interacttran[i]];
    }
    for(int i = 0; i< M; i++)
    {
      for(int j=0; j<N; j++)
      {  
	
	if(a[i][j] < 0 && alreadyused[j] == 0) 
	{
	  lookuptable[i][j1++] = j;
	  alreadyused[j]=1; 
	}
      }
      j1=0; 
    }  
 
 
 
 for(int i = 0; i< M; i++)
   {
    for(int j=0; j<interacttran[i]; j++)
    {
	
	int x = lookuptable[i][j];
	pie[i][j] = conc[i]*prop[x]*k[x];                   // this is the pie vector
	pie_H[i][j] = conc_H[i]*prop[x]*k[x];
	pie_L[i][j] = conc_L[i]*prop[x]*k[x];
	for(int k = 0; k < no_input_spec_per_reac[x]; k++)
	{
	 if ( input_spec_per_reac[x][k] != i )
	 {
	   pie[i][j] = pie[i][j]*conc[input_spec_per_reac[x][k]];
	   pie_H[i][j] = pie_H[i][j]*conc_H[input_spec_per_reac[x][k]];
	   pie_L[i][j] = pie_L[i][j]*conc_L[input_spec_per_reac[x][k]];
	   
	 }
	 if( input_spec_per_reac[x][k] == i && change_input_conc[x][k] == -2)
	 {
	   pie[i][j] = pie[i][j]*(conc[input_spec_per_reac[x][k]]-1)*0.5;
	   pie_H[i][j] = pie_H[i][j]*(conc_H[input_spec_per_reac[x][k]]-1)*0.5;
	   pie_L[i][j] = pie_L[i][j]*(conc_L[input_spec_per_reac[x][k]]-1)*0.5;
	   
	 }
	 if(input_spec_per_reac[x][k] == i && change_input_conc[x][k] == -3)
	 {
	   pie[i][j] = pie[i][j]*(conc[input_spec_per_reac[x][k]]-1)*(conc[input_spec_per_reac[x][k]]-2)/6;
	   pie_H[i][j] = pie_H[i][j]*(conc_H[input_spec_per_reac[x][k]]-1)*(conc_H[input_spec_per_reac[x][k]]-2)/6;
	   pie_L[i][j] = pie_L[i][j]*(conc_L[input_spec_per_reac[x][k]]-1)*(conc_L[input_spec_per_reac[x][k]]-2)/6;
	 }
	}
    }
   }
   
 /*
    for(int i = 0; i< M; i++)
      AND[i] = 0; 
    for(int i = 0; i< M; i++)
   {
     for(int j=0; j<interacttran[i]; j++)
    {
	//cout << "  " << pie[i][j];
	AND[i] = AND[i] + pie[i][j];          // this is the AND (^) vector
	
    }
   }
 */  
/*   for(int i=0; i<M; i++)
     System.out.print(" "+ AND[i]); */

for(int i = 0; i< M; i++)
   {
     for(int j=0; j<interacttran[i]; j++)
    {
	//cout << "  " << pie[i][j];
	sigma[i] = sigma[i] + pie_H[i][j];          // this is the omega vector
	
    }
   }
/*   
   for(int i = 0; i< M; i++)
   {
    sigma[i] = conc_H[i]*AND[i]; 	      // this is the sigma vector
   }
*/   
//   System.out.println(); 
/*   for(int i=0; i<M; i++)
     System.out.print(" "+ sigma[i]); */
   for(int i = 0; i< M; i++)
   {
    totalprop = totalprop + sigma[i];  
   }
 //  System.out.println();
   //System.out.println(totalprop);
  }
  




    
public  void updatedirectedgraph(int selrxn, double[] intervalprop, int[] species_inter)
{
   int l; 
   
  for (int i1 = 0; i1 < no_spec_per_reac[selrxn]; i1++)
  {
     l = spec_per_reac[selrxn][i1]; 
    // System.out.println(" "+conc[l]+" "+conc_L[l]+" "+conc_H[l]);
    
     if((conc[l] != 0 && (conc[l] < conc_L[l] || conc[l] > conc_H[l])) || conc[l] <= 0)
	{ 
     slow_update(l);
     up++;
     totalprop = totalprop - sigma[l]; 
     
       intervalprop[species_inter[l]] = intervalprop[species_inter[l]] - sigma[l];
     
     sigma[l]=0; 
     for(int i=0; i< interacttran[l]; i++)
     {
     int x = lookuptable[l][i];
	pie_L[l][i] = conc_L[l]*k[x];   
	pie_H[l][i] = conc_H[l]*k[x];                  // this is the pie vector
	for(int k = 0; k < no_input_spec_per_reac[x]; k++)
	{
	 if ( input_spec_per_reac[x][k] != l )
	 {
	   pie_H[l][i] = pie_H[l][i]*conc_H[input_spec_per_reac[x][k]];
	   pie_L[l][i] = pie_L[l][i]*conc_L[input_spec_per_reac[x][k]];
	 }
	   
	 if( input_spec_per_reac[x][k] == l && change_input_conc[x][k] == -2)
	 {
	   pie_H[l][i] = pie_H[l][i]*(conc_H[input_spec_per_reac[x][k]]-1)*0.5;
	   pie_L[l][i] = pie_L[l][i]*(conc_L[input_spec_per_reac[x][k]]-1)*0.5;
	   
	 }
	 if(input_spec_per_reac[x][k] == l && change_input_conc[x][k] == -3)
	 {
	   pie_H[l][i] = pie_H[l][i]*(conc_H[input_spec_per_reac[x][k]]-1)*(conc_H[input_spec_per_reac[x][k]]-2)/6;
	   pie_L[l][i] = pie_L[l][i]*(conc_L[input_spec_per_reac[x][k]]-1)*(conc_L[input_spec_per_reac[x][k]]-2)/6;
	 }
	}
	sigma[l] = sigma[l] + pie_H[l][i]; 
     }
     totalprop = totalprop + sigma[l]; 
     
       intervalprop[species_inter[l]] = intervalprop[species_inter[l]] + sigma[l];
    
    int miu = lookuptable[ival][jval];
    int flag = 0;
    for(int k = 0; k < no_input_spec_per_reac[miu]; k++)
    {
		if ( input_spec_per_reac[miu][k] == l )		
		{
			flag = 1;
			break;
		}
	}	
    if (ival != l && flag == 1)
    {
		totalprop = totalprop - pie_H[ival][jval];
		intervalprop[species_inter[ival]] = intervalprop[species_inter[ival]] - pie_H[ival][jval];
		sigma[ival] = sigma[ival] - pie_H[ival][jval];
		pie_H[ival][jval] = conc_H[ival]*conc_H[l]*k[miu];
		totalprop = totalprop + pie_H[ival][jval];
		intervalprop[species_inter[ival]] = intervalprop[species_inter[ival]] + pie_H[ival][jval];
		sigma[ival] = sigma[ival] + pie_H[ival][jval];
	}
	
    }
     
  }
  
} 
    
    
  
  
  
   public static void main(String args[])
   {
     System.out.println(" Enter the values of the number of species and reactions "); 
   BlSSSA obj = new BlSSSA(); 
   
   
  
  System.out.println("Reading File from Java code");
		//Name of the file
		String fileName="xyz.txt", fileName1="k.txt", fileName2="conc.txt";
		try
		{
			//Create object of FileReader
			File fileTm = new File("simulationTime.txt");
			FileReader inputFile = new FileReader(fileName), inputFile1 = new FileReader(fileName1), inputFile2 = new FileReader(fileName2);
			//Instantiate the BufferedReader Class
			BufferedReader bufferReader = new BufferedReader(inputFile), bufferReader1 = new BufferedReader(inputFile1), bufferReader2 = new BufferedReader(inputFile2);
			//Variable to hold the one line data
			String line;
			String[] result;
			
			File file = new File("filenameoutput.csv");
			if (!fileTm.exists()) 
			{
				fileTm.createNewFile();
			}
 			// if file doesnt exists, then create it
			if (!file.exists()) 
			{
				file.createNewFile();
			}
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			FileWriter fwTm = new FileWriter(fileTm.getAbsoluteFile());
			BufferedWriter bwTm = new BufferedWriter(fwTm);
			// Read file line by line and print on the console
			int m11=0, m12=0; 
			while ((line = bufferReader.readLine()) != null)   
			{
				
				//System.out.println(line);
				//result=line.split("\\s");
				result=line.split(",");
				
				for (int i=0; i<obj.N;i++)
				//double value = Double.parseDouble(result[i]);
				{
				int value = Integer.parseInt(result[i]);
				obj.a[m11][m12++]=value;
				//System.out.print("  "+value);
				}
				m11++;
				m12=0; 
				//System.out.println();
			//	bw.write(result[0]);
			//	bw.write("\n");
				
				 
			}
			int count_k=0; 
			while((line = bufferReader1.readLine()) != null)
			{
			//System.out.println(line); 
			obj.k[count_k++] = Double.parseDouble(line);
			}
			//for(int i=0; i<obj.N; i++)
			//System.out.println(obj.k[i]); 
			int[] conc1 = new int[obj.M]; 
			int count_c=0;  
			while((line = bufferReader2.readLine()) != null)
			{
			//System.out.println(line); 
			conc1[count_c++] = Integer.parseInt(line.trim());
			}      
		
  
 
    
    for(int j=0;j<obj.N;j++)
    {
	  for(int i=0;i<obj.M;i++)
	  {
	    if (obj.a[i][j]!= 0)
	    obj.no_spec_per_reac[j]++;  // number of species connected per reaction 
	    if (obj.a[i][j]<0)
	    obj.no_input_spec_per_reac[j]++; // number of input species connected per reaction 
	  }
	    
    }
    for(int i=0; i<obj.N;i++)
    {
      obj.spec_per_reac[i] = new int[obj.no_spec_per_reac[i]];  // species per reaction 
      obj.change_conc[i] = new int[obj.no_spec_per_reac[i]];    // stoichiometry of the species per reaction 
      obj.input_spec_per_reac[i] = new int[obj.no_input_spec_per_reac[i]]; // input species per reaction
      obj.change_input_conc[i] = new int[obj.no_input_spec_per_reac[i]]; // stoichiometry of the input species per reaction 
      
    }
    
    for (int j=0; j< obj.N; j++)
    {	int x=0, x1=0; 
      for ( int i=0; i< obj.M; i++)
      {
	if (obj.a[i][j]!= 0)
	{
	  obj.spec_per_reac[j][x] = i;
	  obj.change_conc[j][x++] = obj.a[i][j]; 
	}
	if (obj.a[i][j] < 0)
	{
	  obj.input_spec_per_reac[j][x1] = i;
	  obj.change_input_conc[j][x1++] = obj.a[i][j]; 
	}
      }
    }
   
     			
    Scanner sc1 = new Scanner(System.in);
  
    System.out.println("Input the fluctuation interval for species populations. For example 0.1 for 10%"); 
  Scanner sc2 = new Scanner(System.in);
  obj.fluc = sc2.nextDouble();
    for (int i=0; i< obj.M; i++)
	{
	    obj.conc[i]=conc1[i]; 
	    obj.slow_update(i);
	    //System.out.println(obj.conc[i]);	
	}

 
    obj.calcpropensity();
    //obj.makeU_3(); 
   // dependencygraph(); 
    System.out.println(" Enter total number of simulation steps");
    int endtime=sc1.nextInt();
    int currenttime=0;
    System.out.println("If you want to write the output in file then, press 1 and if not then, press 0"); 
    int fl = sc1.nextInt();
    long startTime=0, stopTime = 0; long searchtime=0, updatetime=0;
   
    
    
    int intervals = (int) Math.sqrt(obj.M); int[] species_inter = new int[obj.M]; 
   double[] intervalprop = new double[intervals]; int[] interval = new int[intervals]; // intervalprop (BLOCKPROP) is an array that holds the propensity sums of each block 
   System.out.println(" Intervals = "+intervals); 
   int s=0; interval[0]=0;                                         // intervals (NOFBL) is the number of blocks 
   for(int i=0; i<obj.M; i++)					   // interval (BLOCK) is an array which holds the indices of each start species of each block. 
   {								   // species_inter (SPECIES_BLOCK) is an array where ith element holds the index of the block to which the species i belongs 
      //if (i < (s+1)*Math.ceil(obj.M/intervals))
      if (i < (s+1)*obj.M/intervals)
	intervalprop[s] = intervalprop[s] + obj.sigma[i]; 
      else
	{
      	s++; 
      	intervalprop[s] = intervalprop[s] + obj.sigma[i]; 
      	interval[s]=i; 
      	}
      	species_inter[i]=s;
   }
   
   
  
   
   
   
   
   
   
   if(fl == 1)
    {
    System.out.println("Enter the number of simulation steps you want to be written in the file. Steps must be less than total no. of simulation steps "); 
    int inter = sc1.nextInt(), inter1;
    inter = endtime/inter; 
    inter1=inter; int timeplus=0; 
    bw.write("time  ");
    for(int i=0; i< obj.M; i++)
    bw.write("species"+i+"  "); 
    bw.write("\n");
     startTime = System.currentTimeMillis();
    while(currenttime < endtime)
	{
	    //calcpropensity(t, p);
	    long startsearchTime = System.currentTimeMillis();
	 int  selrxn = obj.firingtran(intervalprop, interval, intervals);
	 long stopsearchTime = System.currentTimeMillis();
	 searchtime = searchtime + stopsearchTime - startsearchTime; 
	   
	    obj.transitionfiretime=-Math.log(obj.u)/obj.totalprop;
	    obj.exec(selrxn);
	    long startupdateTime = System.currentTimeMillis();
	    obj.updatedirectedgraph(selrxn, intervalprop, species_inter);
	    long stopupdateTime = System.currentTimeMillis();
	    updatetime = updatetime + stopupdateTime - startupdateTime; 
	    currenttime = currenttime + 1;
	    obj.totaltransitionfiringtime+=obj.transitionfiretime; 
	    if(currenttime == inter1-1)
	    {
	    inter1= inter1+inter; 
	      String strI = String.valueOf(timeplus++);
	      bw.write(strI+"  ");
	      for(int i=0; i<obj.M; i++)
	      {
		 strI = String.valueOf(obj.conc[i]);
		 bw.write(strI+"  ");
	      }
	       bw.write("\n");
	    }
	    //System.out.println(currenttime); 
	}	
	 stopTime = System.currentTimeMillis();
    }
    else
    {
    startTime = System.currentTimeMillis();
    while(currenttime < endtime)
	{
	    //calcpropensity(t, p);
	    long startsearchTime = System.currentTimeMillis();
	 int  selrxn = obj.firingtran(intervalprop, interval, intervals);
	 long stopsearchTime = System.currentTimeMillis();
	 searchtime = searchtime + stopsearchTime - startsearchTime; 
	   // if (selrxn == 0 )
	    //outf << selrxn << " "; 
	    //sparetran++;
	    //obj.transitionfiretime=-Math.log1p(obj.u)/obj.totalprop;
	    obj.transitionfiretime=-Math.log(obj.u)/obj.totalprop;
	    obj.exec(selrxn);
	    long startupdateTime = System.currentTimeMillis();
	    obj.updatedirectedgraph(selrxn, intervalprop, species_inter);
	    long stopupdateTime = System.currentTimeMillis();
	    updatetime = updatetime + stopupdateTime - startupdateTime; 
	    currenttime = currenttime + 1;
	    obj.totaltransitionfiringtime+=obj.transitionfiretime; 
	    //System.out.println(currenttime); 
	}
	stopTime = System.currentTimeMillis();
    }
   
   
   
   
   
   
   
   
   
     
    
   
    System.out.println();
	System.out.println(" "+ obj.conc[0]+" "+ obj.conc[1]+" "+ obj.conc[obj.M-2]+" "+obj.conc[obj.M-1]); 
   
       
      long elapsedTime = stopTime - startTime;
      
      obj.totalprop = 0; obj.deltotalprop=0;  obj.ival=0; obj.jval=0;

      
      System.out.println(" Search time "+searchtime+"ms"); 
      bwTm.write(" Search time "+searchtime+"ms\n"); 
      System.out.println(" Update time "+updatetime+"ms"); 
      bwTm.write(" Update time "+updatetime+"ms\n");
      System.out.println(" Total simulation time "+elapsedTime+"ms");
      bwTm.write(" Total simulation time "+elapsedTime+"ms\n");
      System.out.println("No. of updates "+obj.up);
      bwTm.write("No. of updates "+obj.up+"\n"); 
      System.out.println("No. of trials "+obj.trials);
      bwTm.write("No. of trials "+obj.trials+"\n"); 
      bwTm.close();
      bw.close();
			//Close the buffer reader
			bufferReader.close();
			bufferReader1.close();
			bufferReader2.close();
      }catch(Exception e)
		{
			System.out.println("Error while reading file line by line:"+ e.getMessage());                      
		}
      
    }
    
    
    
   
}
