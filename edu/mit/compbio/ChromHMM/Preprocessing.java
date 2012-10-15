
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

/**
 * This class supports functions to convert either read level data or 
 * signal data to binarized data.
 * The ChromHMM code was written by Jason Ernst 
 */
public class Preprocessing
{

    /**
     * This class takes a gridcontrol array and stores in sumgridcontrol the sum of all values within nflankwidthcontrol bins
     */
    private static void windowSumGrid(int[][][] gridcontrol,int[][][] sumgridcontrol, int nflankwidthcontrol)
    {
	//iterates over chromosome, bin position, and mark
	for (int nchrom = 0; nchrom < gridcontrol.length; nchrom++)
	{
	    int[][] gridcontrol_nchrom = gridcontrol[nchrom];
	    int[][] sumgridcontrol_nchrom = sumgridcontrol[nchrom];

	    for (int nbin = 0; nbin < sumgridcontrol_nchrom.length; nbin++)
	    {
		int[] sumgridcontrol_nchrom_nbin = sumgridcontrol_nchrom[nbin];
		int nstart = Math.max(0,nbin - nflankwidthcontrol);
		int nend = Math.min(nbin + nflankwidthcontrol,gridcontrol_nchrom.length-1);
		for (int nmark = 0; nmark < sumgridcontrol_nchrom_nbin.length; nmark++)
		{
		   int nsum = 0;
                   for (int nrow = nstart; nrow <= nend; nrow++)
		   {
		       int nval = gridcontrol_nchrom[nrow][nmark];
		       if (nval > 0)
		       {
			   nsum += nval;
		       }
		   }
		   sumgridcontrol_nchrom[nbin][nmark] = nsum;
		}
	    }
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts read level data for a cell into integer count information
     * grid is the array for which the read counts will be loaded into
     * bpresent[i] is set to 1 if a there is a read for a chromosome with an index in hmchrom
     * bpresentmark[i] is set to 1 if the mark is present in cell type
     * marks - contains the names of the header marks
     * nshift - is the number of bases a read should be shifted in the 5' to 3' direction of a red
     * nbinsize - is the number of base pairs in a bin
     * bcenterinterval - if true uses the center of the read instead of shifting if the read has already been extended
     * noffsetleft - the amount that should be subtracted from the left coordinate so it is 0-based inclusive
     * noffsetright - the amount that should be subtracted from the right coordinate so it is 0-based inclusive
     * hmfiles - has a mapping of cell\tmark to an acutal read file
     * szcell - the cell we are interested
     * szmarkdir - the directory with bedfiles to read
     * hmchrom - maps chromosome names to an index
     * ninitval - the value that the data should be initialized to
     * szcolfields - a comma delimited string indicating the 0-based columns of the chromosome, start,end,and optionally strand position
     * if null use 0,1,2 for chromosome, start, and end with strand the sixth column or last if fewer
     */
    private static void loadGrid(int[][][] grid,boolean[] bpresent, boolean[] bpresentmarks, String[] marks,int nshift, int nbinsize,
                                 boolean bcenterinterval,int noffsetleft,
				 int noffsetright,HashMap hmfiles, String szcell, String szmarkdir,HashMap hmchrom, 
                                 int ninitval, String szcolfields,boolean bpeaks,boolean bcontrol) throws IOException
    {
	int nummarks = grid[0][0].length;
	//initalizes all values in grid to ninitval
       for (int nchrom = 0; nchrom < grid.length; nchrom++)
       {
          int[][] grid_nchrom = grid[nchrom];
	  int numbins = grid_nchrom.length;
	  for (int nbin = 0; nbin < numbins; nbin++)
	  {
             int[] grid_nchrom_nbin = grid_nchrom[nbin];
	     for (int nmark = 0; nmark < nummarks; nmark++)
	     {
		 grid_nchrom_nbin[nmark] = ninitval;
	     }
          }
       }

       //columns
       int nchromcol=-1;
       int nbegincol=-1;
       int nendcol=-1;
       int nstrandcol= -1;
       int nmaxindex = -1;
       if (szcolfields!=null)
       {
          StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");

          int numtokens = stcolfields.countTokens();
          if ((numtokens <3)||((numtokens<4)&&!bcenterinterval))
          {
             throw new IllegalArgumentException(" invalid number of column fields in "+szcolfields+" expecting 3 or 4 integers");
          }
          nchromcol = Integer.parseInt(stcolfields.nextToken());
          nbegincol = Integer.parseInt(stcolfields.nextToken());
          nendcol = Integer.parseInt(stcolfields.nextToken());
          if (!bcenterinterval)
          {
             nstrandcol = Integer.parseInt(stcolfields.nextToken());
          }

	  nmaxindex = Math.max(nchromcol,Math.max(nbegincol,Math.max(nendcol, nstrandcol)));
       }

       //going through all the mark files in each cell type
       for (int nmark = 0; nmark < nummarks; nmark++)
       {
          boolean bdatafound = false;

	  ArrayList alfiles = (ArrayList) hmfiles.get(szcell+"\t"+marks[nmark]);
	  
	  if ((nummarks == 1)&&(bcontrol))
	  {
	      //this code was added in version 1.04 to handle the situation in which there is a missing mark 
	      //in the first position, but control data can be listed elsewhere
	      //trying to find a listing where control data is available.
	      int nfilemark = 1;
	      while ((alfiles == null)&&(nfilemark < marks.length))
	      {
	         //only one control for the cell looking for a valid one  
	         alfiles = (ArrayList) hmfiles.get(szcell+"\t"+marks[nfilemark]);		
		 nfilemark++;
	      }
	  }

	  if (alfiles == null)
	  {
	     if (bcontrol)
	     {
		 System.out.println("Warning did not find control data for "+szcell+" "+marks[nmark]+" treating as missing");
	     }
	     else
	     {
		 System.out.println("Warning did not find data for "+szcell+" "+marks[nmark]+" treating as missing");
	     }

	     bpresentmarks[nmark] = false;
	     if (!bcontrol)
	     {
		 //slight efficiency improvement here in v1.04
                for (int nchrom = 0; nchrom < grid.length; nchrom++)
                {
                   int[][] grid_nchrom = grid[nchrom];
	           int numbins = grid_nchrom.length;
	           for (int nbin = 0; nbin < numbins; nbin++)
	           {
		       grid_nchrom[nbin][nmark] = -1;		   
		   }
		}
	     } 		               
       	  }
	  else
	  {
	     
	     bpresentmarks[nmark] = true;

             for (int nfile = 0; nfile < alfiles.size(); nfile++)
	     {
		 String szfile = (String) alfiles.get(nfile);
		 BufferedReader brbed = Util.getBufferedReader(szmarkdir+"/"+szfile);
		 String szLine;

	         if (szcolfields != null)
	         {
	            while ((szLine = brbed.readLine())!= null)
	            {
		       String[] szLineA = szLine.split("\\s+");
		       String szchrom = szLineA[nchromcol];

 	               Integer objInt = (Integer) hmchrom.get(szchrom);

		       if (nmaxindex >= szLineA.length)
		       {
			  throw new IllegalArgumentException("Column index "+nmaxindex+" exceeds maximum index "+(szLineA.length-1)+" indicies are 0 based");		           
		       }

		       //if we don't have the chromosome for the read will ignore it
	               if (objInt != null)
	               {
	                  int nchrom = objInt.intValue();
		          int nbin;
		          if (bpeaks)
		          {
			     int nstart = Math.max(0,(Integer.parseInt(szLineA[nbegincol])-noffsetleft)/nbinsize); 
			     int nend = Math.min(grid[nchrom].length-1, (Integer.parseInt(szLineA[nendcol])-noffsetright)/nbinsize);		      

			     for (nbin = nstart; nbin <= nend; nbin++)
	                     {
		                //increment bin count if falls into valid interval
	                        grid[nchrom][nbin][nmark]++;
		                //we do have this chromosome
	                        bpresent[nchrom] = true;			    
			        bdatafound = true;
			     }
			  }
		          else
		          { 
                             if (bcenterinterval)
		             {
		                //uses the center of the interval which is useful if read is already extended
			        nbin = (Integer.parseInt(szLineA[nbegincol])-noffsetleft+Integer.parseInt(szLineA[nendcol])-noffsetright)/(2*nbinsize);
			     }
		             else
		             {
			        String szstrand = szLineA[nstrandcol];
			     	 	      
		                if (szstrand.equals("+"))
		                {		      
			           nbin = (Integer.parseInt(szLineA[nbegincol])-noffsetleft+nshift)/nbinsize; 
                                   //removed one from here may need it for backwards consistency		      		       		      
				}
		                else if (szstrand.equals("-"))
	                        {
		                   nbin = (Integer.parseInt(szLineA[nendcol])-noffsetright-nshift)/nbinsize;		      
				}
		                else
		                {
				    throw new IllegalArgumentException(szstrand+" is an invalid strand!");
				}
			     }
 		   
		             if ((nbin>=0)&&(nbin < grid[nchrom].length))
	                     {
		                //increment bin count if falls into valid interval
	                        grid[nchrom][nbin][nmark]++;
		                //we do have this chromosome
	                        bpresent[nchrom] = true;			    
			        bdatafound = true;
			     }
			  }
		       }
		    }
		    brbed.close();
		 }	     
	         else
	         {
	            while ((szLine = brbed.readLine())!= null)
	            {
	               StringTokenizer st = new StringTokenizer(szLine,"\t ");
		       if (!st.hasMoreTokens())
		       {
		          throw new IllegalArgumentException("Empty line found in "+szmarkdir+"/"+szfile);
		       }
		       String szchrom = st.nextToken();
		       Integer objInt = (Integer) hmchrom.get(szchrom);

		       //if we don't have the chromosome for the read will ignore it
	               if (objInt != null)
	               {
		          int nchrom = objInt.intValue();
 		          if (!st.hasMoreTokens())
		          {
		             throw new IllegalArgumentException("Missing begin coordinate in "+szmarkdir+"/"+szfile);
			  }

			  String szbegin = st.nextToken();

		          if (!st.hasMoreTokens())
		          {
		             throw new IllegalArgumentException("Missing end coordinate in "+szmarkdir+"/"+szfile);
			  }
			  String szend = st.nextToken();
			  int nbin;

		          if (bpeaks)
		          {
			     int nstart = Math.max(0,(Integer.parseInt(szbegin)-noffsetleft)/nbinsize); 
			     int nend = Math.min(grid[nchrom].length-1, (Integer.parseInt(szend)-noffsetright)/nbinsize);		      

			     for (nbin = nstart; nbin <= nend; nbin++)
	                     {
		                //increment bin count if falls into valid interval
	                        grid[nchrom][nbin][nmark]++;
		                //we do have this chromosome
	                        bpresent[nchrom] = true;			    
			        bdatafound = true;
			     }
			  }
		          else
		          {
		             if (bcenterinterval)
		             {
		                //uses the center of the interval which is useful if read is already extended
		                nbin = (Integer.parseInt(szbegin)-noffsetleft+Integer.parseInt(szend)-noffsetright)/(2*nbinsize);
			     }
		             else
		             {
			        if (!st.hasMoreTokens())
			        {
				   throw new IllegalArgumentException("strand column expected, but not found in "+szmarkdir+"/"+szfile);
				}
		                //looks for strand in sixth column or last column if less than six
	                        String szstrand = st.nextToken();
	                        if (st.hasMoreTokens())
	                           szstrand = st.nextToken();
	                        if (st.hasMoreTokens())
	    	                   szstrand = st.nextToken();
			     	 	      
		                if (szstrand.equals("+"))
		                {		      
		                   nbin = (Integer.parseInt(szbegin)-noffsetleft+nshift)/nbinsize; 
                                   //removed one from here may need it for backwards consistency		      		       		      
				}
		                else if (szstrand.equals("-"))
	                        {
		                   nbin = (Integer.parseInt(szend)-noffsetright-nshift)/nbinsize;		      
				}
		                else
		                {
		                   throw new IllegalArgumentException(szstrand+" is an invalid strand!");
				}
			     }
		   
		             if ((nbin>=0)&&(nbin < grid[nchrom].length))
	                     {
		                //increment bin count if falls into valid interval
	                        grid[nchrom][nbin][nmark]++;
		                //we do have this chromosome
	                        bpresent[nchrom] = true;			    
			        bdatafound = true;
			     }
		          }
		       }
		    }	   
		 }
	         brbed.close();
	     }

	     if (!bdatafound)
	     {
	        System.out.println("WARNING not able to load any data for "+szcell+"\t"+marks[nmark]);
	     }
	  }
       }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * This class is used to binarize data both directly and
     * szchromfile is the name of the file with chromsome information
     * it is a two column file with the first being the chromosome and the second being the chromosome length
     * If szcontroldir is null then not using any control
     * szchromlengthfile - two column file which contains the chromosome and strand information 
     * szmarkdir - the directory containing the bed files with the regular mark data
     * szcontroldir - the directory containing the bed files with the control data
     * nflankwidthcontrol - Specifies the number of bins used in both directions to estimate the background; only relevant if control is being used
     * szcellmarkfiletable - tab delimited file each row contains a cell, then mark, then bed file and optionally next a corresponding control bed file
     * nshift - is the number of bases a read should be shifted in the 5' to 3' direction of a red
     * bcenterinterval - if true uses the center of the read instead of shifting if the read has already been extended
     * noffsetleft - the amount that should be subtracted from the left coordinate so it is 0-based inclusive
     * noffsetright - the amount that should be subtracted from the right coordinate so it is 0-based inclusive
     * szoutputsignaldir - if non-null the intermediate signal data will be printed
     * szoutputbinarydir - the directory where the binarized data will be printed
     * szoutputcontroldir - if non-null the intermediate control signal data will be printed
     * dpoissonthresh - the tail probability threshold on the poisson
     * dfoldthresh - the fold threshold required for a present call
     * bcontainsthresh - if true poisson cut off should be highest that still contains dpoissonthresh probability
     * and if false requires strictly greater   
     * npseudocountcontrol -  an integer pseudocount that is uniformay added to every interval to smooth the control data
     * nbinsize - is the number of base pairs in a bin
     * szcolfields - a comma delimited string indicating the 0-based columns of the chromosome, start,end,and optionally strand position
     * if null use 0,1,2 for chromosome, start, and end with strand the sixth column or last if fewer
     */
    public static void makeBinaryDataFromBed(String szchromlengthfile, String szmarkdir, String szcontroldir, int nflankwidthcontrol,String szcellmarkfiletable,
					     int nshift,  boolean bcenterinterval,int noffsetleft, int noffsetright,
                                             String szoutputsignaldir,String szoutputbinarydir, String szoutputcontroldir, 
					     double dpoissonthresh, double dfoldthresh,boolean bcontainsthresh, int npseudocountcontrol,int nbinsize,
					     String szcolfields, boolean bpeaks
                                            ) throws IOException
    {

	//reads in the chromosome length information file
	//the first column of this file is the chromosome and the second is the chromsome length
	BufferedReader brchrom = Util.getBufferedReader(szchromlengthfile);
	String szLine;
	ArrayList allines = new ArrayList();
	while ((szLine = brchrom.readLine())!=null)
	{
	    allines.add(szLine);
	}
	brchrom.close();

	String[] chroms = new String[allines.size()]; //stores the chromosome name
	int[] lengths = new int[chroms.length]; //stores the chromosome length
	HashMap hmchrom = new HashMap(); //stores a mapping of chromomsome to chromosome index
	for (int ni = 0; ni < chroms.length; ni++)
	{
	    StringTokenizer st = new StringTokenizer((String) allines.get(ni),"\t ");
	    if (!st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("empty line found in "+szchromlengthfile);
	    }
	    chroms[ni] = st.nextToken();
	    hmchrom.put(chroms[ni], Integer.valueOf(ni));
	    if (!st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("missing chromosome length for "+allines.get(ni)+" in "+szchromlengthfile);
	    }
	    lengths[ni] = Integer.parseInt(st.nextToken());
	}


	BufferedReader brcellmark = Util.getBufferedReader(szcellmarkfiletable);
	HashSet hscells = new HashSet(); //contains the names of all cell types
	HashSet hsmarks = new HashSet(); //contains the names of all marks
	HashMap hmfiles = new HashMap(); //contains a mapping from (cell type, mark) to the regular data file
	HashMap hmfilescontrol = new HashMap(); //contains a mapping from (cell type, mark) to control file
	HashMap hmfilescellcontrol = new HashMap();//contains a mapping from cell type to a hash set with all its control files

	boolean bcontrol = false; //whether there is control data at all
	String szcontrolfile;

	while ((szLine = brcellmark.readLine())!=null)
	{
	    if (szLine.trim().equals("")) continue;
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
	    if (st.countTokens() < 3)
	    {
		throw new IllegalArgumentException("In "+szcellmarkfiletable+" "+szLine+" has less than 3 columns, expecting at least 3!");
	    }
	    String szcell = st.nextToken();
	    String szmark = st.nextToken();
	    String szfile = st.nextToken();
	    if (st.hasMoreTokens())
	    {
		//we have control data
		szcontrolfile = st.nextToken();
		bcontrol = true;

		//was hmfiles in version 1.00
	        ArrayList alfilescontrol = (ArrayList) hmfilescontrol.get(szcell+"\t"+szmark);
	        if (alfilescontrol == null)
	        {
		   alfilescontrol = new ArrayList();
		   hmfilescontrol.put(szcell+"\t"+szmark,alfilescontrol);
		}
		alfilescontrol.add(szcontrolfile);

		HashSet hscellcontrol = (HashSet) hmfilescellcontrol.get(szcell);
		if (hscellcontrol == null)
		{
		    hscellcontrol = new HashSet();
		    hmfilescellcontrol.put(szcell, hscellcontrol);
		}		
		
	        hscellcontrol.add(szcontrolfile);		
	    }
	    hscells.add(szcell);
	    hsmarks.add(szmark);

	    ArrayList alfiles = (ArrayList) hmfiles.get(szcell+"\t"+szmark);
	    if (alfiles == null)
	    {
		alfiles = new ArrayList();
	        hmfiles.put(szcell+"\t"+szmark,alfiles);
	    }
	    alfiles.add(szfile);
	}
	brcellmark.close();
	
	//loads all the marks in hsmarks into the array marks and then sorts it
	int nummarks = hsmarks.size();
	int numcontrolmarks=-1;
	String[] marks = new String[nummarks];
	int nmarkindex = 0;
	Iterator itrmarks = hsmarks.iterator();
	while (itrmarks.hasNext())
	{
	    marks[nmarkindex] = (String) itrmarks.next();
	    nmarkindex++;
	}
	Arrays.sort(marks);	


	int[][][] grid = new int[chroms.length][][]; //first dimension is chromosome, second is number of full size bins, third is mark
	int[][][] gridcontrol = null; //similiar to grid but only created if 
	int[][][] sumgridcontrol = null; //a smoothed version of grid control
	boolean[] bpresentcontrol = null; //flags which chromosome we actually have control data for
	boolean[] bpresentmarkscontrol = null;

	if (bcontrol)
	{
	    //control data is indicated going to allocate memory for it
           gridcontrol = new int[chroms.length][][];
           sumgridcontrol = new int[chroms.length][][];
	   bpresentcontrol = new boolean[chroms.length];
	   bpresentmarkscontrol = new boolean[nummarks];
	}

	for (int ni = 0; ni < chroms.length; ni++)
	{
	    //allocating the full memory for the real data
	    grid[ni] = new int[lengths[ni]/nbinsize][nummarks];
	}

	boolean[] bpresent = new boolean[chroms.length]; //chromosome we actually have read data for
	boolean[] bpresentmarks = new boolean[nummarks];

	Iterator itrcells = hscells.iterator();
	while (itrcells.hasNext())
	{
	    //going through each declared cell type
	    String szcell = (String) itrcells.next();
	    HashSet hscellcontrol = (HashSet) hmfilescellcontrol.get(szcell);

	    boolean bcontrolfile;
            if (hscellcontrol == null)
	    {
		//no control data for this cell type
		bcontrolfile = false;
	    }
	    else
	    {
		bcontrolfile = true;
		if (hscellcontrol.size()==1)
		{
		    //we have one control for all marks
		    numcontrolmarks = 1;
		}
		else
		{
		    //will allocate the full memory for all marks
		    numcontrolmarks = nummarks;
		}
	    }


	    //loading data for the cell type
	    loadGrid(grid,bpresent,bpresentmarks,marks,nshift,nbinsize,bcenterinterval,noffsetleft,
		     noffsetright,hmfiles,szcell,szmarkdir,hmchrom,0,szcolfields,bpeaks,false);
	    if (bcontrolfile)
	    {
	       if ((gridcontrol[0] == null)||(gridcontrol[0][0].length !=numcontrolmarks))
	       {
		   //reallocate if changing array size
		   //allowed to go between single and matched 
	          for (int ni = 0; ni < chroms.length; ni++)
	          {
	             gridcontrol[ni] = new int[lengths[ni]/nbinsize][numcontrolmarks];
	             sumgridcontrol[ni] = new int[lengths[ni]/nbinsize][numcontrolmarks];
	          }
	       }

	       //we have control data loading cell type data for that
	       loadGrid(gridcontrol,bpresentcontrol,bpresentmarkscontrol,marks,nshift,nbinsize,bcenterinterval,noffsetleft,noffsetright,
                        hmfilescontrol,szcell,szcontroldir,hmchrom,npseudocountcontrol,szcolfields,bpeaks,true);
	    }
	    
	    if (szoutputsignaldir!=null)
	    {
		//printing signal is requested
	       for (int nchrom = 0; nchrom < chroms.length; nchrom++)
	       {
		  if (bpresent[nchrom])
		  {
		      //we have signal for this chromosome
		      String szfile = szoutputsignaldir+"/"+szcell+"_"+chroms[nchrom]+"_signal.txt";
		      System.out.println("Writing to file "+szfile);
		      PrintWriter pw = new PrintWriter(szfile);
		      pw.println(szcell+"\t"+chroms[nchrom]);
		      //outputs mark header
                      for (int nmark = 0; nmark < marks.length-1; nmark++)
                      {
		         pw.print(marks[nmark]+"\t");
		      }
		      pw.println(marks[marks.length-1]);

		      //outputs mark signal data
		      int[][] grid_nchrom = grid[nchrom];
 	              for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
                      { 
			 int[] grid_nchrom_nbin = grid_nchrom[nbin];
	                 for (int nmark = 0; nmark < grid_nchrom_nbin.length-1; nmark++)
	                 {
	                    pw.print(grid_nchrom_nbin[nmark]+"\t");
		         }
		         pw.println(grid_nchrom_nbin[marks.length-1]);
		      }
		      pw.close();	       
		  }
	       }


	       if ((bcontrolfile)&&(szoutputcontroldir!=null))
	       {
                  for (int nchrom = 0; nchrom < chroms.length; nchrom++)
                  {
	             if (bpresent[nchrom])
	             {
	                //we have signal for this chromosome
			String szfile = szoutputcontroldir+"/"+szcell+"_"+chroms[nchrom]+"_controlsignal.txt";
		        System.out.println("Writing to file "+szfile);
	                PrintWriter pw = new PrintWriter(szfile);
			pw.println(szcell+"\t"+chroms[nchrom]);
	                //outputs mark header
			if (numcontrolmarks == 1)
			{
			    pw.println("Control");
			}
			else
			{
                           for (int nmark = 0; nmark < marks.length-1; nmark++)
                           {
       	                      pw.print(marks[nmark]+"\t");
	                   }
	       	           pw.println(marks[marks.length-1]);
			}

		        //outputs mark signal data
		        int[][] gridcontrol_nchrom = gridcontrol[nchrom];
 	                for (int nbin = 0; nbin < gridcontrol[nchrom].length; nbin++)
                        {  
	        	   int[] gridcontrol_nchrom_nbin = gridcontrol_nchrom[nbin];
			   if (gridcontrol_nchrom_nbin.length == 1)
			   {
		              pw.println(gridcontrol_nchrom_nbin[0]-npseudocountcontrol);
			   }
			   else
			   {
	                      for (int nmark = 0; nmark < gridcontrol_nchrom_nbin.length-1; nmark++)
	                      {
				  pw.print((gridcontrol_nchrom_nbin[nmark]-npseudocountcontrol)+"\t");
		              }
		              pw.println(gridcontrol_nchrom_nbin[marks.length-1]-npseudocountcontrol);
			   }
		        }
		        pw.close();	       
		     }
		  }		  
	       }
	    }

	    int nummarks_m1 = nummarks - 1;

	    if (bcontrolfile)
	    {
	       //binarization will be based on control data

	       //smoothing control data
	       windowSumGrid(gridcontrol,sumgridcontrol,nflankwidthcontrol);

	       //determiming thresholds for each mark and background depth
	       int[][] thresholds = null;

	       if (!bpeaks)
	       {
                  thresholds = determineMarkThresholdsFromBinnedDataArrayAgainstControl(grid,sumgridcontrol,
				        bpresent,bpresentcontrol,dpoissonthresh,dfoldthresh,bcontainsthresh);
	       }

	       for (int nchrom = 0; nchrom < chroms.length; nchrom++)
               {
       	          if ((bpresent[nchrom])&&(bpresentcontrol[nchrom]))
	          {
		      String szfile = szoutputbinarydir+"/"+szcell+"_"+chroms[nchrom]+"_binary.txt";
		      System.out.println("Writing to file "+szfile);
		      PrintWriter pw = new PrintWriter(szfile);
	       	     //we have both primary and control data for the mark
		      pw.println(szcell+"\t"+chroms[nchrom]);

                     for  (int nmark = 0; nmark < nummarks_m1; nmark++)
		     {
      	                pw.print(marks[nmark]+"\t");
	             }
	             pw.println(marks[nummarks_m1]);
		     
		     int[][] grid_nchrom = grid[nchrom];
	      	     int[][] sumgridcontrol_nchrom = sumgridcontrol[nchrom];
 	             for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
                     {  
	                int[] grid_nchrom_nbin = grid_nchrom[nbin];
			int[] sumgrid_nchrom_nbin = sumgridcontrol_nchrom[nbin];

	                for (int nmark = 0; nmark < nummarks_m1; nmark++)
	                {
			    int ncontrolval;
			    if (numcontrolmarks == 1)
			    {
				ncontrolval = sumgrid_nchrom_nbin[0];
			    }
			    else
			    {
				ncontrolval = sumgrid_nchrom_nbin[nmark];
			    }
	      		   //printing one if count exceeds background threshold

			   if (!bpresentmarks[nmark])
			   {
			       pw.print("2\t");
			   }
		           else if ((bpeaks&&(grid_nchrom_nbin[nmark]>=1)) ||(!bpeaks&&(thresholds[nmark][ncontrolval] <= grid_nchrom_nbin[nmark])))
		           {
	                       pw.print("1\t");
	      	           }
		           else
		           {
	                       pw.print("0\t");
		           }
		        }

			int ncontrolval;
		        if (numcontrolmarks == 1)
		        {
			   ncontrolval = sumgrid_nchrom_nbin[0];
		        }
		        else
		        {
		       	   ncontrolval = sumgrid_nchrom_nbin[nummarks_m1];
		        }


			if (!bpresentmarks[nummarks_m1])
		        {
	                   pw.println("2");
			}
		        else if ((bpeaks&&(grid_nchrom_nbin[nummarks_m1]>=1))||(!bpeaks&&(thresholds[nummarks_m1][ncontrolval] <= grid_nchrom_nbin[nummarks_m1])))
	                {
	                   pw.println("1");
	      	        }
		        else
	                {
	      	           pw.println("0");
		        }
		     }
		     pw.close();
		  }
	          else
	          {
	             if (bpresent[nchrom]&&!bpresentcontrol[nchrom])
		     {
	                System.out.println("WARNING for "+szcell+" "+chroms[nchrom]+" regular data is found but not control data");
		     }
		     else if (!bpresent[nchrom]&&bpresentcontrol[nchrom])
		     {
		        System.out.println("WARNING for "+szcell+" "+chroms[nchrom]+" control data is found but not regular data");
		     } 
		  }	       
	       }    
	    }
	    else
	    {
              //going to use background thresholds based on uniform background model
	       double[] thresholds = null;

	       if (!bpeaks)
	       {
                  thresholds = determineMarkThresholdsFromBinnedDataArray(grid,bpresent,dpoissonthresh,dfoldthresh,bcontainsthresh);
	       }
		
	      for (int nchrom = 0; nchrom < chroms.length; nchrom++)
	      {
	         if (bpresent[nchrom])
	         {
		    String szfile = szoutputbinarydir+"/"+szcell+"_"+chroms[nchrom]+"_binary.txt";
		    System.out.println("Writing to file "+szfile);
	            PrintWriter pw = new PrintWriter(szfile);
		    pw.println(szcell+"\t"+chroms[nchrom]);
                    for (int nmark = 0; nmark < marks.length-1; nmark++)
                    {
	               pw.print(marks[nmark]+"\t");
	            }
	            pw.println(marks[nummarks_m1]);
		
	            int[][] grid_nchrom = grid[nchrom];
                    for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
                    {  
	               //printing 1 if signal has met data threshold and 0 otherwise
	              int[] grid_nchrom_nbin = grid_nchrom[nbin];
                      for (int nmark = 0; nmark < nummarks_m1; nmark++)
                      {
			 if (!bpresentmarks[nmark])
		         {
	                    pw.print("2\t");
			 }
                         else if ((bpeaks&&(grid_nchrom_nbin[nmark]>=1))||(!bpeaks&&(thresholds[nmark] <= grid_nchrom_nbin[nmark])))
	       	         {
	                    pw.print("1\t");
	        	 }
		         else
		         {
	                    pw.print("0\t");
	                 }
		      }

		      if (!bpresentmarks[nummarks_m1])
		      {
	                 pw.println("2");
		      }
		      else if ((bpeaks&&(grid_nchrom_nbin[nummarks_m1]>=1))||(!bpeaks&&(thresholds[nummarks_m1] <= grid_nchrom_nbin[nummarks_m1])))
	              {
                         pw.println("1");
	       	      }
	              else
	              {
	                 pw.println("0");
		      }
		   }
	           pw.close();
		 }	       
	      }	       
	   }	   
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * grid - the integer data values from which to determine the poisson cutoffs
     * gridcontrol - the control data which the thresholds will be relative to
     * bpresent - a vector indicating which indicies of grid to include in the analysis
     * bpresentcontrol - a vector indicating which indicies of gridcontrol to include in the analysis
     * dpoissonthresh - the tail probability threshold on the poisson
     * dfoldthresh - the fold threshold required for a present call
     * bcontainsthresh - if true poisson cut off should be highest that still contains dpoissonthresh probability
     * and if false requires strictly greater     
     **/
    public static int[][] determineMarkThresholdsFromBinnedDataArrayAgainstControl(int[][][] grid, int[][][] gridcontrol,
										   boolean[] bpresent, boolean[] bpresentcontrol,
										   double dpoissonthresh, double dfoldthresh,boolean bcontainsthresh)
    {
     
       double dcumthreshold = 1-dpoissonthresh;

       int nummarks= grid[0][0].length;
       int numcontrolmarks = gridcontrol[0][0].length;

       //stores the the total number of reads for each mark and its matched control
       long[] sumtags = new long[nummarks];
       long[] sumtagscontrol = new long[nummarks];

       //stores the thresholds for each mark and background value
       int[][] thresholds = new int[nummarks][];

       //stores the maximum control value found for each mark
       int[] maxcontrol = new int[nummarks];

       //stores which background control values have been found for each mark
       HashSet[] hscontrol = new HashSet[nummarks];
       for (int nmark = 0; nmark < nummarks; nmark++)
       {
	   hscontrol[nmark] = new HashSet();
       }


       for (int nchrom = 0; nchrom < grid.length; nchrom++)
       {
	   if ((bpresent[nchrom])&&(bpresentcontrol[nchrom]))
	   {
	      int[][] grid_nchrom = grid[nchrom];
	      int[][] gridcontrol_nchrom = gridcontrol[nchrom];

	      for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
	      {
	         int[] grid_nchrom_nbin = grid[nchrom][nbin];
	         int[] gridcontrol_nchrom_nbin = gridcontrol[nchrom][nbin];
                 for (int nmark = 0; nmark < nummarks; nmark++)
                 {
		    int nval = grid_nchrom_nbin[nmark];
		    int ncontrolval;

		    if (numcontrolmarks == 1)
		    {
                       ncontrolval = gridcontrol_nchrom_nbin[0];
		    }
		    else
		    {
                       ncontrolval = gridcontrol_nchrom_nbin[nmark];
		    }

		    if (ncontrolval > maxcontrol[nmark])
	            {
	       	       maxcontrol[nmark] = ncontrolval;
		    }
		    hscontrol[nmark].add(Integer.valueOf(ncontrolval));
	            sumtags[nmark] += grid_nchrom_nbin[nmark];
	            sumtagscontrol[nmark] += ncontrolval;		            
		 }
	      }
	   }
       }	  
      
        for (int nmark = 0; nmark < sumtags.length; nmark++)
        {
	    //computing threshold for each mark
	   thresholds[nmark] = new int[maxcontrol[nmark]+1];
	   int[] thresholds_nmark = thresholds[nmark];

           //determine the relative enrichment for real reads versus the local expected
	   double davgratio = sumtags[nmark]/(double)sumtagscontrol[nmark];

	   //sets a background of 0 threshold to 1
	   thresholds_nmark[0] = 1;
	   
	   //going through each background value
           for (int nbackground = 1; nbackground < maxcontrol[nmark]; nbackground++)
           {
	       if(hscontrol[nmark].contains(Integer.valueOf(nbackground)))
	       {
		   //only compute the background threshold for values we observed

		   //expected number of reads is the local background number of reads times the global
		   //readdepth enrichment for sumtags
                  double dlambda = davgratio*nbackground;

	          double dcum = 0;
	          short nthresh = 0;
	          double dlogfactorial = 0;
 
	          while (dcum <= dcumthreshold)
                  {
		     double dprob = Math.exp(Math.log(dlambda)*nthresh-dlambda- dlogfactorial);
	             dcum += dprob;
	             nthresh++;
	 	     dlogfactorial += Math.log(nthresh);
	          }

		  if (bcontainsthresh)
		  {
		      //decreasing to include the dpoissonthreshold porbability
                     nthresh--;
		  }

                  thresholds_nmark[nbackground] = Math.max(Math.max((int) Math.ceil(dfoldthresh*dlambda),nthresh),1);

		  if (ChromHMM.BVERBOSE)
		  {
                     System.out.println("Modification\t"+nmark+"\t"+nbackground+"\t"+maxcontrol[nmark]+"\t"+thresholds_nmark[nbackground]+"\t"+dlambda);
		  }
	       }
	   }
	}
	return thresholds;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * grid - the integer data values from which to determine the poisson cutoffs
     * bpresent - a vector indicating which indicies of grid to include in the analysis
     * dpoissonthresh - the tail probability threshold on the poisson
     * dfoldthresh - the fold threshold required for a present call
     * bcontainsthresh - if true poisson cut off should be highest that still contains dpoissonthresh probability
     * and if false requires strictly greater
     **/
    public static double[] determineMarkThresholdsFromBinnedDataArray(int[][][] grid, boolean[] bpresent, 
								      double dpoissonthresh,double dfoldthresh,boolean bcontainsthresh)
    {
       double dcumthreshold = 1-dpoissonthresh;

       double[] sumtags =null;
       double[] thresholds = null;
       int nummarks= grid[0][0].length;
       String szcell;

       int ntotallocs = 0;

       sumtags = new double[nummarks];
       thresholds = new double[nummarks];

       //computes for those chromosomes considered present
       //the total number of locations in ntotallocs
       //and the total signal for the mark sumtags[nmark]
       for (int nchrom = 0; nchrom < grid.length; nchrom++)
       {
	   if (bpresent[nchrom])
	   {
	      int[][] grid_nchrom = grid[nchrom];
	      for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
	      {
	         int[] grid_nchrom_nbin = grid[nchrom][nbin];
                 for (int nmark = 0; nmark < nummarks; nmark++)
                 {		    
	            sumtags[nmark] += grid_nchrom_nbin[nmark];
		 }	      
	      }
	      ntotallocs += grid_nchrom.length;
	   }
        }	   

       for (int nj = 0; nj < sumtags.length; nj++)
       {
	   //computes expected number of reads in a bin
          double dlambda = sumtags[nj]/(double) ntotallocs;
	  double dcum = 0;
	  short nthresh = 0;

	  double dlogfactorial = 0;
	  while (dcum <= dcumthreshold)
          {
	      //lambda^x times e^(-lambda) / (x!)
	      //taking the exponential of the log term
	      double dprob = Math.exp(Math.log(dlambda)*nthresh-dlambda-dlogfactorial);
	      dcum += dprob;
	      nthresh++;
	      dlogfactorial += Math.log(nthresh);
	  }

	  if (bcontainsthresh)
	  {
	      //want to ensure  dpoissonthresh is included with the cutoff
             nthresh--;
	  }

	  //threshold set to greater of poisson based and fold threshold
          thresholds[nj] = Math.max(Math.max((int) Math.ceil(dfoldthresh*dlambda),nthresh),1);
	  if (ChromHMM.BVERBOSE)
	  {
             System.out.println("Threshold\t"+nj+"\t"+thresholds[nj]+"\t"+sumtags[nj]+"\t"+((int) Math.ceil(dfoldthresh*dlambda))+"\t"+ntotallocs);
	  }
       }
       return thresholds;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Binarize from signal data with a control
     * Control data can either be one control data for all marks or a matched control for each mark
     * szbinneddataDIR - directory containing the files with the signal data. files of the same cell type have the same prefix before the '_'
     * szcontrolDIR - directory containing the files with the control data. these files should have the same name as their matched file 
     * in szbinneddataDIR
     * szoutputDIR - directory to which the binarized data will be written. The '_binary.txt' will be appended
     * dpoissonthresh - the tail probability threshold on the poisson
     * dfoldthresh - the fold threshold required for a present call
     * bcontainsthresh - if true poisson cut off should be highest that still contains dpoissonthresh probability
     * and if false requires strictly greater               
     * nflankwidthcontrol -  Specifies the number of bins used in both directions to estimate the background; only relevant if control is being used
     * npseudocountcontrol - an integer pseudocount that is uniformay added to every interval to smooth the control data
     */
    public static void makeBinaryDataFromSignalAgainstControl(String szbinneddataDIR, String szcontrolDIR, String szoutputDIR,
							      double dpoissonthresh,double dfoldthresh, boolean bcontainsthresh, int nflankwidthcontrol, 
                                                              int npseudocountcontrol) throws IOException
    {
       int nummarks=-1;
       int nummarkscontrol=-1;
       File dir = new File(szbinneddataDIR);
       String[] allfiles = dir.list();
       int nfilecount = 0;
       for (int nfile = 0; nfile < allfiles.length; nfile++)
       {
	   if (allfiles[nfile].contains("_signal"))
	   {
	       nfilecount++;
	   }
       }
       String[] signalchromfiles = new String[nfilecount];
       int nfileindex = 0;
       for (int nfile = 0; nfile < allfiles.length; nfile++)
       {
	   if (allfiles[nfile].contains("_signal"))
	   {
	       signalchromfiles[nfileindex] = allfiles[nfile];
	       nfileindex++;
	   }
       }

       //maps cell type chrom to a control file
       HashMap hmcontrol = new HashMap();

       File dircontrol = new File(szcontrolDIR);
       String[] allfilescontrol = dircontrol.list();
       for (int nfile = 0; nfile < allfilescontrol.length; nfile++)
       {
	   if (allfilescontrol[nfile].contains("_controlsignal"))
	   {
	       BufferedReader brcontrol = Util.getBufferedReader(szcontrolDIR+"/"+allfilescontrol[nfile]);
	       String szheader = brcontrol.readLine();
	       if (szheader == null)
	       {
		   throw new IllegalArgumentException(szcontrolDIR+"/"+allfilescontrol[nfile]+" is empty!");
	       }
	       StringTokenizer st = new StringTokenizer(szheader,"\t");
	       if (st.countTokens() != 2)
	       {
		   throw new IllegalArgumentException(szcontrolDIR+"/"+allfilescontrol[nfile]+" header line must have two columns delimited by a tab found only one in "+
						      szheader);
	       }
	       String szkey = st.nextToken()+"\t"+st.nextToken();
	       hmcontrol.put(szkey,allfilescontrol[nfile]);
	       brcontrol.close();
	   }
       }

       HashMap hmcellsToIndex = new HashMap();
       String szcell;


       //goes through all signal files
       //storing for each cell type, an ArrayList containing the indicies of files for that cell type
       for (int nfile = 0; nfile < signalchromfiles.length; nfile++)
       {
	   BufferedReader br = Util.getBufferedReader(szbinneddataDIR+"/"+signalchromfiles[nfile]);
	   String szLine = br.readLine();
	   String szcurrcell;
	   if (szLine == null)
	   {
	       throw new IllegalArgumentException(signalchromfiles[nfile]+" is empty!");
	   }
	   StringTokenizer st = new StringTokenizer(szLine,"\t"); 
	   szcurrcell = st.nextToken();
	   br.close();

          ArrayList al = (ArrayList) hmcellsToIndex.get(szcurrcell);

          if (al == null)
	  {
	      al = new ArrayList();
	      hmcellsToIndex.put(szcurrcell, al);
	  }
          al.add(Integer.valueOf(nfile));
       }


       Iterator itrcells = hmcellsToIndex.entrySet().iterator();

       while (itrcells.hasNext())
       {
	   Map.Entry pairs = (Map.Entry) itrcells.next();
	   //going through all the cells
	   szcell = (String) pairs.getKey();
	   ArrayList al = (ArrayList) pairs.getValue();

	   int numchrompercell = al.size();
	   int ntotallocs = 0;

	   int[][][] grid = new int[numchrompercell][][];
	   int[][][] gridcontrol = new int[numchrompercell][][];
	   int[][][] sumgridcontrol = new int[numchrompercell][][];
	   String[] chroms = new String[numchrompercell];
	   String szHeaderLine2= null;
	   for (int nchrom = 0; nchrom < numchrompercell; nchrom++)
	   {
	       //going through all the files for the cell type

	       String szfilename = signalchromfiles[((Integer) al.get(nchrom)).intValue()];
	       BufferedReader br = Util.getBufferedReader(szbinneddataDIR+"/"+szfilename);

	       String szLine;
	       String szHeaderLine1 = br.readLine();
	       szHeaderLine2 = br.readLine();

	       if (szHeaderLine1 == null)
	       {
		   throw new IllegalArgumentException(szcontrolDIR+"/"+szfilename+" is empty!");
	       }
	       StringTokenizer stheader = new StringTokenizer(szHeaderLine1,"\t");
	       stheader.nextToken();
	       if (!stheader.hasMoreTokens())
	       {
		   throw new IllegalArgumentException("Only found one entry for line "+szHeaderLine1+" in file "+szbinneddataDIR+"/"+szfilename
						      +" expecting 2");
	       }
	       String szchrom = stheader.nextToken();
	       chroms[nchrom] = szchrom;
	       String szcontrolfilename = (String) hmcontrol.get(szcell+"\t"+szchrom);
	       if (szcontrolfilename == null)
	       {
		   throw new IllegalArgumentException("No control data found for "+szcell+"\t"+szchrom);
	       }
	       BufferedReader brcontrol = Util.getBufferedReader(szcontrolDIR+"/"+szcontrolfilename);
	       String szHeaderLineControl1 = brcontrol.readLine();
	       String szHeaderLineControl2 = brcontrol.readLine();

	       StringTokenizer st = new StringTokenizer(szHeaderLine2,"\t");
	       int ncurrnummarks = st.countTokens();

	       if (szHeaderLineControl2 == null)
	       {
		   throw new IllegalArgumentException(szcontrolDIR+"/"+szcontrolfilename+" is missing lines!");
	       }
	       StringTokenizer stcontrol = new StringTokenizer(szHeaderLineControl2,"\t");
	       int ncurrnummarkscontrol = stcontrol.countTokens();
	       if ((ncurrnummarkscontrol != 1) && (ncurrnummarkscontrol != ncurrnummarks))
	       {
		   throw new IllegalArgumentException("Number of marks in control file "+szcontrolDIR+"/"+szfilename+
                                       " is "+ncurrnummarkscontrol+" expecting "+ncurrnummarks+" or 1");
	       }

	       if (nchrom == 0)
	       {
		   nummarks = ncurrnummarks;
		   nummarkscontrol = ncurrnummarkscontrol;
	       }
	       else
	       {
                  if (ncurrnummarks != nummarks)
	          {
		     throw new IllegalArgumentException("Number of marks ("+nummarks+ ") in "+szbinneddataDIR+"/"+
                                                       szfilename+" is inconsistent with previously found # of "+nummarks);
		  }
	          else if (ncurrnummarkscontrol != nummarkscontrol)
	          {
		     throw new IllegalArgumentException("Number of marks ("+nummarkscontrol+ ") in "+szcontrolDIR+"/"+szfilename+
                                                      " is inconsistent with previously found # of "+nummarks);
		  }
	       }


	       //storing the read sum for all marks
	       int nchromlocs = 0;
 	       while ((szLine = br.readLine())!=null)
	       {
	          nchromlocs++;
	       }
	       br.close();
	       ntotallocs += nchromlocs;

	       //allocationg memory for the data
	       grid[nchrom] = new int[nchromlocs][nummarks];
	       gridcontrol[nchrom] = new int[nchromlocs][nummarkscontrol];
	       sumgridcontrol[nchrom] = new int[nchromlocs][nummarkscontrol];

	       int nbin = 0;

	       br = Util.getBufferedReader(szbinneddataDIR+"/"+szfilename);
	       br.readLine(); //gets rid of the headers
	       br.readLine();

	       int[][] grid_nchrom = grid[nchrom];
	       int[][] gridcontrol_nchrom = gridcontrol[nchrom];
 	       while ((szLine = br.readLine())!=null)
	       {
		  int[] grid_nchrom_nbin = grid_nchrom[nbin];
		  int[] gridcontrol_nchrom_nbin = gridcontrol_nchrom[nbin];
	          st = new StringTokenizer(szLine,"\t");
	          for (int nmark = 0; nmark < nummarks; nmark++)
                  {
		      //reading in the regular data
		      grid_nchrom_nbin[nmark] = Integer.parseInt(st.nextToken());		      
		  }

		  String szLineControl = brcontrol.readLine();
		  if (szLineControl == null)
		  {
		      throw new IllegalArgumentException("The number of lines in the control file "+
			     szcontrolDIR+"/"+szfilename+" does not match that in the signal file "+szbinneddataDIR+"/"+szcontrolfilename); 
		  }
	          st = new StringTokenizer(szLineControl,"\t");
	          for (int nmark = 0; nmark < nummarkscontrol; nmark++)
                  {
		      //reading in the control data
		      gridcontrol_nchrom_nbin[nmark] = Integer.parseInt(st.nextToken())+npseudocountcontrol;		      
		  }
		  nbin++;
	       }
	       br.close();
	       brcontrol.close();
	   }

           //binarization will be based on control data
           //smoothing control data
           windowSumGrid(gridcontrol,sumgridcontrol,nflankwidthcontrol);

	   boolean[] bpresent = new boolean[grid.length];
	   boolean[] bpresentcontrol = new boolean[sumgridcontrol.length];

	   for (int ni = 0; ni < bpresent.length; ni++)
           {
	      //all chromosomes are valid
	      bpresent[ni] = true;
	      bpresentcontrol[ni] = true;
	   }

           //determiming thresholds for each mark and background depth
           int[][] thresholds = determineMarkThresholdsFromBinnedDataArrayAgainstControl(grid,sumgridcontrol,
											     bpresent,bpresentcontrol,dpoissonthresh,dfoldthresh,bcontainsthresh);

	   for (int nchrom = 0; nchrom < grid.length; nchrom++)
           {
	      String szchrom = chroms[nchrom];
	      String szfile = szoutputDIR+"/"+szcell+"_"+szchrom+"_binary.txt";
	      System.out.println("Writing to file "+szfile);
              PrintWriter pw = new PrintWriter(new FileWriter(szfile));	 
	      
	      int nummarks_m1 = nummarks - 1;
	      pw.println(szcell+"\t"+szchrom);
	      pw.println(szHeaderLine2);

	      int[][] grid_nchrom = grid[nchrom];
              int[][] sumgridcontrol_nchrom = sumgridcontrol[nchrom];
 	      for (int nbin = 0; nbin < grid_nchrom.length; nbin++)
              {   
	         int[] grid_nchrom_nbin = grid_nchrom[nbin];
	         int[] sumgridcontrol_nchrom_nbin = sumgridcontrol_nchrom[nbin];

	         for (int nmark = 0; nmark < nummarks_m1; nmark++)
	         {
		    int ncontrolval;
		    if (nummarkscontrol == 1)
	            {
		       //always use the first control mark column if only one column
		       ncontrolval = sumgridcontrol_nchrom_nbin[0];
		    }
		    else
	            {
		       ncontrolval = sumgridcontrol_nchrom_nbin[nmark];
		    }

	      	    //printing one if count exceeds background determined threshold
		    if (grid_nchrom_nbin[nmark] == -1)
		    {
			pw.print("2\t");
		    }
                    else if (thresholds[nmark][ncontrolval] <= grid_nchrom_nbin[nmark])
	            {
                       pw.print("1\t");
            	    }
	            else
	            {
                       pw.print("0\t");
   	            }
		 }

	         int ncontrolval;
	         if (nummarkscontrol == 1)
	         {
	            ncontrolval = sumgridcontrol_nchrom_nbin[0];
	         }
	         else
	         {
	            ncontrolval = sumgridcontrol_nchrom_nbin[nummarks_m1];
	         }


		 if (grid_nchrom_nbin[nummarks_m1] == -1)
	         {
	       	    pw.println("2");
	         }
                 else if (thresholds[nummarks_m1][ncontrolval] <= grid_nchrom_nbin[nummarks_m1]) 
	         {
	            pw.println("1");
	         }
	         else
	         {
                    pw.println("0");
                 }
	      }
	      pw.close();
           }
       }       
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Binarizes data from signal data using a uniform background model
     * szbinneddataDIR - the directory containing the binned signal data
     * szoutputDIR - the directory to which the binarized data should be written. Files have the same file name as in the
     * szbinneddatDIR with '_binary.txt' extension added
     * dpoissonthresh - the tail probability threshold on the poisson
     * dfoldthresh - the fold threshold required for a present call
     * bcontainsthresh - if true poisson cut off should be highest that still contains dpoissonthresh probability
     * and if false requires strictly greater     
     **/
    public static void makeBinaryDataFromSignalUniform(String szbinneddataDIR, String szoutputDIR,
                                                       double dpoissonthresh,double dfoldthresh, boolean bcontainsthresh) throws IOException
    {
	//this computes the binarization without storing in main memory all the data

       double dcumthreshold = 1-dpoissonthresh;
       double[] sumtags =null;
       int[] thresholds = null;

       int nummarks=-1;
       File dir = new File(szbinneddataDIR);

       String[] allfiles = dir.list();
       int nfilecount = 0;
       for (int nfile = 0; nfile < allfiles.length; nfile++)
       {
	   if (allfiles[nfile].contains("_signal"))
	   {
	       nfilecount++;
	   }
       }
       String[] signalchromfiles = new String[nfilecount];
       int nfileindex = 0;
       for (int nfile = 0; nfile < allfiles.length; nfile++)
       {
	   if (allfiles[nfile].contains("_signal"))
	   {
	       signalchromfiles[nfileindex] = allfiles[nfile];
	       nfileindex++;
	   }
       }
       HashMap hmcellsToIndex = new HashMap();
       String szcell;

       //goes through all signal files
       //storing for each cell type, an ArrayList containing the indicies of files for that cell type
       for (int nfile = 0; nfile < signalchromfiles.length; nfile++)
       {
	   BufferedReader br = Util.getBufferedReader(szbinneddataDIR+"/"+signalchromfiles[nfile]);
	   String szLine = br.readLine();
	   if (szLine == null)
	   {
	       throw new IllegalArgumentException(szbinneddataDIR+"/"+signalchromfiles[nfile]+" does not contain any data!");
	   }
	   StringTokenizer st = new StringTokenizer(szLine,"\t"); 
	   String szcurrcell = st.nextToken();
	   br.close();

           ArrayList al = (ArrayList) hmcellsToIndex.get(szcurrcell);

           if (al == null)
	   {
	      al = new ArrayList();
	      hmcellsToIndex.put(szcurrcell, al);
	   }
           al.add(Integer.valueOf(nfile));
       }

       Iterator itrcells = hmcellsToIndex.entrySet().iterator();

       while (itrcells.hasNext())
       {
	   Map.Entry pairs = (Map.Entry) itrcells.next();
	   //going through all the cells
	   szcell = (String) pairs.getKey();
	   ArrayList al = (ArrayList) pairs.getValue();
	   int numchrompercell = al.size();
	   int ntotallocs = 0;
	   sumtags = null;

	   for (int nindex = 0; nindex < numchrompercell; nindex++)
	   {
	       //going through all the files for the cell type

	       String szfilename = signalchromfiles[((Integer) al.get(nindex)).intValue()];
	       BufferedReader br = Util.getBufferedReader(szbinneddataDIR+"/"+szfilename);
	       br.readLine();//get rid of header
	       String szLine = br.readLine();
	       if (szLine == null)
	       {
		   throw new IllegalArgumentException(szbinneddataDIR+"/"+szfilename+" is empty!!!");
	       }
	       StringTokenizer st = new StringTokenizer(szLine,"\t");
	       nummarks = st.countTokens();

 	       if (sumtags == null)
               {
		   //found out the number of marks allocating space
		   sumtags = new double[nummarks];
		   thresholds = new int[nummarks];
	       }
	       else if (nummarks != sumtags.length)
	       {
		   throw new IllegalArgumentException("Number of marks ("+nummarks+ ") in "+szfilename+" is inconsistent with previously found # of"+sumtags.length);
	       }

	       //storing the read sum for all marks
 	       while ((szLine = br.readLine())!=null)
	       {
	          ntotallocs++;
	          st = new StringTokenizer(szLine,"\t");
	          for (int nj = 0; nj < nummarks; nj++)
                  {
		      double dval =  Double.parseDouble(st.nextToken());		      
		      sumtags[nj] += dval;
		  }
	       }
	       br.close();
	   }	   
       
           for (int nj = 0; nj < sumtags.length; nj++)
           {
	       //determing the cutoff for each mark
	       double dlambda = sumtags[nj]/(double) ntotallocs;

	       double dcum = 0;
	       short nthresh = 0;
	       double dlogfactorial = 0; 

	       while (dcum <= dcumthreshold)
	       {
		  double dprob = Math.exp(Math.log(dlambda)*nthresh-dlambda-dlogfactorial);
	          dcum += dprob;	
	          nthresh++;
		  dlogfactorial += Math.log(nthresh);
	       }

	       if (bcontainsthresh)
	       {
	          //want to ensure  dpoissonthresh is included with the cutoff
                  nthresh--;
	       }

               thresholds[nj] = Math.max(Math.max((int) Math.ceil(dfoldthresh*dlambda),nthresh),1);
	       if (ChromHMM.BVERBOSE)
	       {
                  System.out.println("Modification\t"+szcell+"\t"+nj+"\t"+thresholds[nj]+"\t"+sumtags[nj]+
				     "\t"+((int) Math.ceil(dfoldthresh*dlambda))+"\t"+ntotallocs);
	       }
	   }


	   for (int nindex = 0; nindex < numchrompercell; nindex++)
	   {
	       //re-reads the input file and outputs the data based on the binarization threshold
	       String szfilename = signalchromfiles[((Integer) al.get(nindex)).intValue()];
	       BufferedReader br = Util.getBufferedReader(szbinneddataDIR+"/"+szfilename);
	       String szChromCellLine = br.readLine();
	       String szMarkLine = br.readLine();
	       if (szChromCellLine == null)
	       {
		   throw new IllegalArgumentException(szbinneddataDIR+"/"+szfilename+" is empty!");
	       }
	       StringTokenizer st = new StringTokenizer(szChromCellLine,"\t");
	       if (st.countTokens() != 2)
	       {
		   throw new IllegalArgumentException(szbinneddataDIR+"/"+szfilename+" header line must have a two columns delimited by a tab found "+
						      szChromCellLine);
	       }
	       String szfile = szoutputDIR+"/"+st.nextToken()+"_"+st.nextToken()+"_binary.txt";
	       System.out.println("Writing to file "+szfile);
	       PrintWriter pw = new PrintWriter(new FileWriter(szfile));
	       pw.println(szChromCellLine);
	       pw.println(szMarkLine);
	       String szLine;
	       while ((szLine = br.readLine())!=null)
	       {
	          st = new StringTokenizer(szLine,"\t");
	          for (int ncol = 0; ncol < nummarks-1; ncol++)
	          {
		     double dval = Double.parseDouble(st.nextToken());
		     if (dval == -1)
	             {
			 pw.print("2\t");
		     }
		     else if (thresholds[ncol] <= dval)
		     {
		        pw.print("1\t");			       
		     }
		     else
		     {
		        pw.print("0\t");
		     }
		  }

		  double dval = Double.parseDouble(st.nextToken());
		  if (dval == -1)
		  {
		     pw.println("2");
		  }
		  else if (thresholds[nummarks-1] <= dval)
		  {
		     pw.println("1");
	          }
		  else
	          {
	      	     pw.println("0");
		  }
	       }
	       br.close();
	       pw.close();
	   }
       }
    }
}
