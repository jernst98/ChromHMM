/**
 * ChromHMM - automating chromatin state discovery and characterization 
 * Copyright (C) 2008-2012 Massachusetts Institute of Technology
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 **/
package edu.mit.compbio.ChromHMM;

import java.io.*;
import java.util.*;

public class NestedEliminateInitialize
{

    /**
     * Implements a procedure to generate a nested initialization based on state pruning
     * The ChromHMM code was written by Jason Ernst 
     */

    /**
     * This is used to generate a set of nested models by pruning states from a high scoring model. This can be used as one heuristic to initialize model parameters 
     * to get a roughly comparable set of models across different number of states while biasing the learning procedure to avoid redundant or non-representative states. 
     * States are greedily pruned from the highest scoring model. The criteria to a prune a state is that it has the least impact on the total distance to the nearest 
     * remaining state for all the states in other models. If there is only one model, then distance to other states in that model are considered. By default euclidean 
     * distance for emission parameters is used, but there is an option to use the total correlation instead.
     * szinputdir is the name of the directory containing model files to eliminate. Only files with prefix 'model_' in the directory are used. 
     * szoutput dir is the name of the output directory where the eliminated model files will go
     * beuclid is true if euclidean distance should be used
     */
    public static void nestedEliminateInitialize(String szinputdir, String szoutputdir,boolean beuclid) throws IOException
    {
        int numcols = 0;
        int nummodels = 0;

        File inputdir = new File(szinputdir);
        String[] files = inputdir.list();

        double dbestlikelihood = Double.NEGATIVE_INFINITY;
	String szmainfile="";

	int nbestnumstates=-1;
	int nbestmodel=-1;
	char chbestorder='0';

	HashMap hmNameToID = new HashMap();
	boolean bfound = false;

        for (int nfile = 0; nfile < files.length; nfile++)
	{
	    //goes through all the files counting how many begin with model_
	    //also for the first found determines the number of marks
	    //and a mapping from marks to IDs

	   if (files[nfile].startsWith("model_"))
	   {
	       if (!bfound)
	       {
		   BufferedReader brfile = Util.getBufferedReader(szinputdir+"/"+files[nfile]);
		   String szLine;
		   numcols = 0;

		   while ((szLine = brfile.readLine())!= null)
		   {
		      if (szLine.startsWith("emission"))
		      {
			  StringTokenizer stemiss = new StringTokenizer(szLine,"\t");
			  stemiss.nextToken();
			  stemiss.nextToken();
			  stemiss.nextToken();
			  String szmark = stemiss.nextToken();
			  
			  Object objInt = hmNameToID.get(szmark);
			  if (objInt ==null)
			  {
			      hmNameToID.put(szmark, Integer.valueOf(numcols));			  
			      numcols++;
			  }
		       }
		   }
		   brfile.close();
		   bfound = true;
	       }
	       nummodels++;
	   }
	}

	if (nummodels < 1)
	{
	    throw new IllegalArgumentException("No models found in directory "+szinputdir);
        }

        double[][][] modelemissions = new double[nummodels][][];
        int nmodel = 0;

        for (int nfile = 0; nfile < files.length; nfile++)
	{
           //going through all the model files and loading in the emission parameters
	    //also recording the model with the best score

	   if (files[nfile].startsWith("model_"))
	   {
	       //assuming file
	       String szfile = files[nfile];

	       BufferedReader br =  Util.getBufferedReader(szinputdir+"/"+szfile);

	       String szheader = br.readLine();
	       if (szheader == null)
		   throw new IllegalArgumentException(szinputdir+"/"+szfile+" is empty!");
	       StringTokenizer stheader = new StringTokenizer(szheader,"\t ");
	       int numstates = Integer.parseInt(stheader.nextToken());
	       stheader.nextToken();
	       char chorder = stheader.nextToken().charAt(0); //order type
	       double dcurrlikelihood = Double.parseDouble(stheader.nextToken());

	       String szLine;
	       StringTokenizer st;

	       modelemissions[nmodel] = new double[numstates][numcols];
	       int ncount = 0;
	       while ((szLine = br.readLine())!=null)
	       {

	          if (szLine.startsWith("emissionprobs"))
		  {
		      st = new StringTokenizer(szLine,"\t");
		      st.nextToken();
		      int nstate = Integer.parseInt(st.nextToken())-1;
		      st.nextToken();
		      String szmark = st.nextToken();
		      Integer intMark = (Integer) hmNameToID.get(szmark);
		      int nmark = -1;
		      if (intMark == null)
		      {
			  throw new IllegalArgumentException(szmark+" not found in "+files[nfile]);
		      }
		      else
		      {
			  nmark = ((Integer) intMark).intValue();
		      }
		      String szbucket = st.nextToken();
		      if (szbucket.equals("1"))
		      {
		         modelemissions[nmodel][nstate][nmark] = Double.parseDouble(st.nextToken());
			 ncount++;
		      }
		  }
	       }
	       br.close();

	       if (ncount != (numstates*numcols))
	       {
		   throw new IllegalArgumentException(files[nfile]+" had "+ncount+" emission values while expecting "+(numstates*numcols));
	       }

	       if (dcurrlikelihood > dbestlikelihood)
	       {
		   dbestlikelihood = dcurrlikelihood;	       
		   nbestmodel = nmodel;
		   nbestnumstates = numstates;
		   chbestorder = chorder;
		   szmainfile = files[nfile];
	       }

	       nmodel++;
	   }
	}
	   
        //storing the parameters for the best model
	double[] bestprobinit = new double[nbestnumstates];
        double[][] besttransitionprobs = new double[nbestnumstates][nbestnumstates];
	ArrayList[] emissionsline = new ArrayList[nbestnumstates];
	   
	BufferedReader brfile =  Util.getBufferedReader(szinputdir+"/"+szmainfile);
        brfile.readLine();
	for (int ni = 0; ni < bestprobinit.length; ni++)
        {
	   String szinitline = brfile.readLine();
	   if (szinitline == null)
	   {
	       throw new IllegalArgumentException(szinputdir+"/"+szmainfile+" is missing lines!");
	   }
	   StringTokenizer st = new StringTokenizer(szinitline,"\t");
           st.nextToken();
           st.nextToken();
           bestprobinit[ni] = Double.parseDouble(st.nextToken());
        }

	for (int ni = 0; ni < nbestnumstates; ni++)
        {
	   for (int nj = 0; nj < nbestnumstates; nj++)
           {
	      String sztransitionline = brfile.readLine();
	      if (sztransitionline == null)
	      {
		  throw new IllegalArgumentException(szinputdir+"/"+szmainfile+" is missing lines!");
	      }
	      StringTokenizer st = new StringTokenizer(sztransitionline,"\t");
	      st.nextToken();
	      st.nextToken();
	      st.nextToken();
	      besttransitionprobs[ni][nj] = Double.parseDouble(st.nextToken());
	   }
	}

        for (int ni = 0; ni < nbestnumstates; ni++)
        {
           emissionsline[ni] = new ArrayList();
	}

	String szLine;
        while ((szLine = brfile.readLine())!=null)
	{
	   StringTokenizer st = new StringTokenizer(szLine,"\t");
           st.nextToken();
	   int nj = Integer.parseInt(st.nextToken())-1;
           emissionsline[nj].add(st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken());
	}
        brfile.close();


        double dmaxcorrouterouter = -Integer.MAX_VALUE;
        int nstatesleft = nbestnumstates;
        boolean[] elim = new boolean[nbestnumstates];
	double[][] modelemissions_nbestmodel = modelemissions[nbestmodel];

        while (nstatesleft > 2)
	{
	    //while there are more than 2 states in the model

	   int nmaxelimstate = -1;
	   double dmaxtotsum = -Double.MAX_VALUE;
	   for (int nelimstate = 0; nelimstate < elim.length; nelimstate++)
	   {
	       //will consider eliminating each state to see effect on total representation
	      if (!elim[nelimstate])
	      {
		  elim[nelimstate] = true;
		  double dtotsum =0;
		  for (int nc = 0; nc < modelemissions.length; nc++)
		  {
		      //will iterate through each sampled models

		     if ((nbestmodel != nc)||(modelemissions.length==1))
		     {
		        //going through each state in the other sampled modes
		        for (int nd = 0; nd < modelemissions[nc].length; nd++)
			{
			    double dmaxagreeval =-Integer.MAX_VALUE;
			    double[] modelemissions_ncnd = modelemissions[nc][nd];
			    //comparing to each state in the best model not eliminated
			    for (int nb = 0; nb < modelemissions_nbestmodel.length; nb++)
			    {			       
			       if (!elim[nb])
			       {
			          //state has not been elminiated
				   double dagreeval;
				   if (beuclid)
				   {
                                      dagreeval = -Util.euclid(modelemissions_nbestmodel[nb],modelemissions_ncnd);
				   }
				   else
				   {
                                      dagreeval = Util.correlation(modelemissions_nbestmodel[nb],modelemissions_ncnd);
				   }

				  if (dagreeval > dmaxagreeval)
			          {
				      //this is the closest state found so far
				     dmaxagreeval = dagreeval;
				  }
			       }
			    }
			    //add to the total distance for this emission vector the closest distance
			    dtotsum += dmaxagreeval;
			}
		     }
		  }

		  //the elimination of this state was temporary
		  elim[nelimstate] = false;
		  if (dtotsum > dmaxtotsum)
	          {
		      dmaxtotsum = dtotsum;
		      nmaxelimstate = nelimstate;
		  }
	      }
	   }

	   if (ChromHMM.BVERBOSE)
	   {
	      System.out.println("****\t"+nstatesleft+"\t"+nmaxelimstate+"\t"+dmaxtotsum);
	   }
	   elim[nmaxelimstate] = true;
	   nstatesleft--;
	   //nmaxelim state is the state to eliminate

	   PrintWriter pw = new PrintWriter(szoutputdir+"/elim_"+nstatesleft+"_"+szmainfile);
	   String szLinein;

	   int nelim = nstatesleft;
	   pw.println(nelim+"\t"+numcols+"\t"+chbestorder);

	   //readjusting the inital probabilities
	   //uniformly assign the initial probabilities of eliminated states
	   int nj = 0;
	   double dextra = 0;
	   double[] probinit = new double[nelim];
	   for (int ni = 0; ni < elim.length; ni++)
           {
	       if (!elim[ni])
	       {
		   probinit[nj] = bestprobinit[ni];
		   nj++;
	       }
	       else
	       {
		   dextra +=  bestprobinit[ni];
	       }
	   }

	   dextra= dextra/nelim;
	   nj = 0;
	   for (int ni = 0; ni < elim.length; ni++)
	   {
	      if (!elim[ni])
	      {
		  pw.println("probinit\t"+nj+"\t"+(probinit[nj]+dextra));
		  nj++;
	      }
	   }

	   //readjusting the transition probabilities
	   int nk = 0;
	   for (int na = 0; na < elim.length; na++)
	   {
	      if (!elim[na])
	      {
		  nj = 0;
		  dextra = 0;
		  double[] transitionprob = new double[nelim];
		  for (int ni = 0; ni < elim.length; ni++)
		  {
		      if (!elim[ni])
		      {
			  transitionprob[nj] = besttransitionprobs[na][ni];
			  nj++;
		      }
		      else
		      {
			  dextra +=  besttransitionprobs[na][ni];
		      }
		  }

		  dextra= dextra/nelim;
		  nj = 0;
		  for (int ni = 0; ni < elim.length; ni++)
		  {
		     if (!elim[ni])
		     {
			 pw.println("transitionprobs\t"+nk+"\t"+nj+"\t"+(transitionprob[nj]+dextra));
			 nj++;
		     }
		  }
		  nk++;
	      }
	   }

	   //outputing emission parameters for those states which have not been eliminated
	   nk = 0;
	   for (int na = 0; na < elim.length; na++)
           {
	      if (!elim[na])
	      {
		 int numlines = emissionsline[na].size();
	         for (int ni = 0; ni < numlines; ni++)
		 {
		     pw.println("emissionprobs\t"+nk+"\t"+emissionsline[na].get(ni));	  
		 }
		 nk++;
	      }
	   }
	   pw.close();
	}
    }
}
