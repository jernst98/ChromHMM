
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
import java.text.*;
import java.awt.*;
import org.tc33.jheatchart.HeatChart;


/**
 * This class supports post analysis for the ChromHMM segmentations
 * The ChromHMM code was written by Jason Ernst 
 */
public class StateAnalysis
{



    /**
     * A record for storing a segmentation
     */
    static class SegmentRec
    {
	String szchrom; // the chromosome
	int nbegin; //coordinate of start of the segmentation
	int nend; //coordinate of the end of the segmentation
	short slabel; //the segmentationlabel

	SegmentRec(String szchrom, int nbegin, int nend,short slabel)
	{
	    this.szchrom = szchrom;
	    this.nbegin =nbegin;
	    this.nend = nend;
	    this.slabel = slabel;
	}
    }

    ////////////////////////////////////////////////////////////////////////////////

    /**
     * Stores all the emission parameters for a model
     * corresponding to szfilename and also explicitly the number of states
     */
    static class RecEmissionFile
    {
	String szfilename;
	int numstates;
	double[][] emissionparams;
    }



    /////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Compares two emission files and orders them so the one with fewer states comes first
     * and if there is a tie then by the file name
     */
    public static class RecEmissionFileCompare implements Comparator, Serializable
    {
	public int compare(Object o1, Object o2)
	{
	    RecEmissionFile r1 = (RecEmissionFile) o1;
	    RecEmissionFile r2 = (RecEmissionFile) o2;
	    if (r1.numstates < r2.numstates)
	    {
		return -1;
	    }
	    else if (r1.numstates > r2.numstates)
	    {
		return 1;
	    }
	    else
	    {
		return (r1.szfilename.compareTo(r2.szfilename));
	    }
	}
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static class Interval
    {
	int nstart;
	int nend;

	Interval(int nstart, int nend)
	{
	    this.nstart = nstart;
	    this.nend = nend;
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static class IntervalCompare implements Comparator, Serializable
    {
	public int compare(Object o1, Object o2)
        {
	    Interval i1 = (Interval) o1;
	    Interval i2 = (Interval) o2;

	    if (i1.nstart < i2.nstart)
	    {
		return -1;
	    }
	    else if (i1.nstart > i2.nstart)
	    {
		return 1;
	    }
	    else if (i1.nend < i2.nend)
	    {
		return -1;
	    }
	    else if (i1.nend > i2.nend)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////



    static void makeEnrichmentHeatMap(double[][] heatmapfold, String[] xlabels, String[] ylabels,String szoutfileprefix,Color theColor,
				      String sztitle,String szxaxis,String szyaxis) throws IOException
     {
       HeatChart map = new HeatChart(heatmapfold);

       map.setTitle(sztitle);
       map.setXAxisLabel(szxaxis);
       map.setYAxisLabel(szyaxis);

       map.setChartMargin(100);
       map.setXValues(xlabels);
       map.setYValues(ylabels);
       map.setHighValueColour(theColor);
       map.setAxisValuesFont(new Font("SansSerif",0,20));
       map.setAxisLabelsFont(new Font("SansSerif",0,22));
       map.setTitleFont(new Font("SansSerif",0,24));
       map.saveToFile(new File(szoutfileprefix+".png"));	
       Util.printImageToSVG(map, szoutfileprefix+".svg");
       System.out.println("Writing to file "+szoutfileprefix+".png");
       System.out.println("Writing to file "+szoutfileprefix+".svg");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     /**
      *
      * Computes an enrichment based on a set of posterior files
      * szposteriordir - the directory with the posterior files to compute enrichments for
      * szcell - the cell type for which to compute the enrichment
      * szinputcoorddir - the directory containing the coordinates for enrichment
      * szinputcoordlist - if non-null only computes enrichments for file names in the specified
      * noffsetleft - the amount that should be subtracted from the left coordinate so it is 0-based inclusive
      * noffsetright - the amount that should be subtracted from the right coordinate so it is 0-based inclusive
      * nbinsize -  the binsize used in the segmentation
      * bcenter - if true only uses the center of the interval for determining enrichments
      * bunique - only counts a bin once, this is not applicable if bbaseres is true or busesignal is true
      * busesignal - if true then uses signal information about the interval if available
      * szcolfields - comma delimited list of the 0-based columns of chromosome,start,end, and if busesignal is true the signal column
      * bbaseres - if true gives enrichment at the base resolution opposed to the bin resolution
      * szoutfile - the file to which to output the 
      * bcolscaleheat - if true scales the heatmap individual for each column otherwise uses the same for all columns
      * theColor - theColor to use for the heatmap

      * default is chromsome is field 1, start is field 2, end is field 3, and amount is field 4
      * if szcolfields is non-null reads chrom, start, end, signal if busesignal field is true
      *
      */
      public static void enrichmentPosterior(String szposteriordir,String szcell,String szinputcoorddir, String szinputcoordlist,
					    int noffsetleft, int noffsetright, int nbinsize,
					    boolean bcenter,boolean bunique, boolean busesignal,String szcolfields,boolean bbaseres, 
					     String szoutfile,boolean bcolscaleheat, Color theColor,String sztitle,String szlabelmapping) throws IOException
    {

	String szLine;
	char chorder='0';

	if (busesignal)
	{
	    bunique = false;
	}

	File posteriordir = new File(szposteriordir);
	String[] posteriorfiles = posteriordir.list();

	int numposteriorstates = 0;

	int nfirstindex=0;
	boolean bfirst = true;


	String[] files;

	if (szinputcoordlist == null)
        {
	   File dir = new File(szinputcoorddir);

	   if (dir.isDirectory())
	   {
              files = dir.list();
              Arrays.sort(files);
	      szinputcoorddir += "/";
	   }
	   else
	   {
	      files = new String[1];
	      files[0] = szinputcoorddir;
	      szinputcoorddir = "";
	   }
	}
	else
	{
	    szinputcoorddir += "/";
	    //loads in the input coords list
	    BufferedReader brfiles = Util.getBufferedReader(szinputcoordlist);

	    ArrayList alfiles = new ArrayList();
	    while ((szLine = brfiles.readLine())!=null)
	    {
		alfiles.add(szLine);
	    }
	    brfiles.close(); 
	    files = new String[alfiles.size()];
	    for (int nfile = 0; nfile < files.length; nfile++)
	    {
		files[nfile] = (String) alfiles.get(nfile);
	    }
	}


	double dsumlabel = 0;
	double[] tallylabel = null;
	double[][] tallyoverlaplabel = null;
        double[] dsumoverlaplabel = new double[files.length];

	boolean bposteriorfound = false;

	for (int nfile = 0; nfile < posteriorfiles.length; nfile++)
	{
	    String szposteriorfiles_nfile = posteriorfiles[nfile];

	    //going through the posterior files
	    if (szposteriorfiles_nfile.contains("_posterior"))
	    {
		
		BufferedReader brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		String szHeader = brposterior.readLine();
		if (szHeader == null)
		    throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" is empty!");
		StringTokenizer st =new StringTokenizer(szHeader,"\t");
		String szcurrcell = st.nextToken();
		if ((!szcurrcell.equals(szcell))&&(!szcell.equals("")))
		{
		    brposterior.close();
		}
		else
		{
		   bposteriorfound = true;
		   String szchrom = st.nextToken();

		   int numlines = 0;
		   szLine = brposterior.readLine();
		   if (szLine == null)
		   {
		       throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" is missing lines!");
		   }
	           st = new StringTokenizer(szLine,"\t");
	           int numcurrstates = st.countTokens();
	           if (bfirst)
	           {
		      chorder = st.nextToken().charAt(0);
		      bfirst = false;
		      nfirstindex = nfile;
		      numposteriorstates = numcurrstates;
		      tallylabel = new double[numposteriorstates];
		      tallyoverlaplabel = new double[files.length][numposteriorstates];
		   }
		   else if (numposteriorstates != numcurrstates)
		   {
	              throw new IllegalArgumentException("Number of states "+numcurrstates+" in "+szposteriorfiles_nfile+
							" does not match number of states "+numposteriorstates+" in "+posteriorfiles[nfirstindex]);
		   }

	           while ((szLine = brposterior.readLine())!=null)
	           {
		      numlines++;
		   }
		   brposterior.close();	       
	    	
		   //now know the number of states and lines for the posterior file
	           float[][] posterior = new float[numlines][numposteriorstates];
		
		   //loading in the posterior data and tallying how frequently each state has occured
	           brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		   brposterior.readLine();
		   brposterior.readLine();
	           int nline = 0;
	           while ((szLine = brposterior.readLine())!=null)
	           {
                      st = new StringTokenizer(szLine,"\t");
		      float[] posterior_nline = posterior[nline];
                      for (int nstate = 0; nstate < numposteriorstates; nstate++)
                      {
		         float fval = Float.parseFloat(st.nextToken());
		         posterior_nline[nstate] = fval;
		         tallylabel[nstate] += fval;
		      }
		      nline++;
		   }
		   brposterior.close();
 
		   for (int ncoordfile = 0; ncoordfile < files.length; ncoordfile++)
		   {		
		      double[] tallyoverlaplabel_ncoordfile = tallyoverlaplabel[ncoordfile];

		      int nchromindex = 0;
		      int nstartindex = 1;
		      int nendindex = 2;
		      int nsignalindex = 3;

	              if (szcolfields  != null)
	              {
		         //gets the start and end coordinates
		         StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
		         nchromindex = Integer.parseInt(stcolfields.nextToken());
		         nstartindex = Integer.parseInt(stcolfields.nextToken());
		         nendindex = Integer.parseInt(stcolfields.nextToken());

	                 if (busesignal)
	                 {
		            nsignalindex = Integer.parseInt(stcolfields.nextToken());
			 }
		      }

		      BufferedReader brcoords = Util.getBufferedReader(szinputcoorddir+files[ncoordfile]);
		      if (bunique)
		      {
			 ArrayList alrecs = new ArrayList();
	                 while ((szLine = brcoords.readLine())!=null)
	                 {
		            if (szLine.trim().equals("")) continue;
		            String[] szLineA = szLine.split("\\s+");
	                    String szreadchrom = szLineA[nchromindex];
		            if (szreadchrom.equals(szchrom))
		            {
	                       String szcurrchrom = szLineA[nchromindex];		             
		               int nbeginactual =Integer.parseInt(szLineA[nstartindex])-noffsetleft;
			       int nendactual =Integer.parseInt(szLineA[nendindex])-noffsetright;
			       if (bcenter)
		               {
			           nbeginactual = (nbeginactual+nendactual)/2;
			           nendactual = nbeginactual;
			       }
			       alrecs.add(new Interval(nbeginactual,nendactual));			       
			    }
			 }

		         Object[] alrecA = alrecs.toArray();
		         Arrays.sort(alrecA,new IntervalCompare());

			 boolean bclosed = true;
		         int nintervalstart = -1;
		         int nintervalend = -1;
		         boolean bdone = false;

		         for (int nindex = 0; ((nindex <= alrecA.length)&&(alrecA.length>0)); nindex++)
		         {
			    int ncurrstart=-1;
			    int ncurrend=-1;

		            if (nindex == alrecA.length)
		            {
			       bdone = true;
			    }
		            else 
		            {
		               ncurrstart = ((Interval) alrecA[nindex]).nstart;
		               ncurrend = ((Interval) alrecA[nindex]).nend;
                               if (nindex == 0)
		               {
			          nintervalstart = ncurrstart;
			          nintervalend = ncurrend;
			       }
		               else if (ncurrstart <= nintervalend)
		               {
			          //this read is still in the active interval
		                  //extending the current active interval interval
		                 if (ncurrend > nintervalend)
		                 {
		                    nintervalend = ncurrend;
			         }
			       }		    
		               else 
		               {
		                  //just finished the current active interval
			          bdone = true;
			       }
			    }

		            if (bdone)
		            {		      						
			       int nbegin = nintervalstart/nbinsize;
			       int nend = nintervalend/nbinsize;
			       if (nend >= posterior.length)
			       {
		                  nend = posterior.length - 1;
			       }

			       for (int nbin = nbegin; nbin <= nend; nbin++)
			       {
			          float[] posterior_nbin = posterior[nbin];
			          for (int nstate = 0; nstate < numposteriorstates; nstate++)
			          {
			             tallyoverlaplabel_ncoordfile[nstate] += posterior_nbin[nstate];
			          }
			       }

		               if (bbaseres)
		               {  
		                  //dbeginfrac represents the fraction of bases the nbegin interval
		                  //which came after the actual nbeginactual
	                          double dbeginfrac = (nintervalstart - nbegin*nbinsize)/(double) nbinsize;
			   
			          //dendfrac represents the fraction of bases after the end position in the interval
			          double dendfrac = ((nend+1)*nbinsize-nintervalend-1)/(double) nbinsize;

			          if ((nbegin < posterior.length)&&(dbeginfrac>0))
		                  { 
				     //only removing if it would of made it into the actual posterior counts
				     float[] posterior_nbegin = posterior[nbegin];
			             for (int nstate = 0; nstate < numposteriorstates; nstate++)
			             {
		                        tallyoverlaplabel_ncoordfile[nstate]-=dbeginfrac*posterior_nbegin[nstate];
			      	     }
				  }

                                  if ((nend < posterior.length)&&(dendfrac >0))
		                  {
				     //only removing if it would of made it into the actual posterior counts		          
				     float[] posterior_nend = posterior[nend];
			             for (int nstate = 0; nstate < numposteriorstates; nstate++)
			             {
			                tallyoverlaplabel_ncoordfile[nstate]-=dendfrac*posterior_nend[nstate];
				     }				      
				  }		   
				  nintervalstart = ncurrstart;
				  nintervalend = ncurrend;
			          bdone = false;				  		  
			       }
			    }
			 }		      
		      }
		      else
		      {
	                 while ((szLine = brcoords.readLine())!=null)
	                 {
		            if (szLine.trim().equals("")) continue;
		            String[] szLineA = szLine.split("\\s+");
	                    String szreadchrom = szLineA[nchromindex];
		            if (szreadchrom.equals(szchrom))
		            {
			       int nbeginactual =Integer.parseInt(szLineA[nstartindex])-noffsetleft;
			       int nbegin = nbeginactual/nbinsize;

			       int nendactual =Integer.parseInt(szLineA[nendindex])-noffsetright;
			       int nend = nendactual/nbinsize;

			       double damount;
	                       if ((busesignal)&&(nsignalindex < szLineA.length))
	                       {
	       	                  damount = Double.parseDouble(szLineA[nsignalindex]);
			       }
	                       else
	                       {
		                  damount = 1;
			       }

		               if (bcenter)
	                       {
			          int ncenter = (nbeginactual+nendactual)/(2*nbinsize);
	                          if (ncenter < posterior.length)
	                          {
			             float[] posterior_ncenter = posterior[ncenter];
		                     
				     //allowing duplicates
			             for (int nstate = 0; nstate < numposteriorstates; nstate++)
			             {  
		                        tallyoverlaplabel_ncoordfile[nstate]+=damount*posterior_ncenter[nstate];
			             }				     
				  }
			       }			    			    		      
	                       else
	                       {
			          //allowing duplicate counts
			          if (nend >= posterior.length)
			          {
			             nend = posterior.length - 1;
				  }

	                          for (int nindex = nbegin; nindex <= nend; nindex++)
	                          {
			             float[] posterior_nindex = posterior[nindex];
			             for (int nstate = 0; nstate < numposteriorstates; nstate++)
			             {
		  	                tallyoverlaplabel_ncoordfile[nstate]+=damount*posterior_nindex[nstate];
				     }
				  }

			          if (bbaseres)
			          { 
				     //the fractions of the begin and end intervals not actually included
	                             double dbeginfrac = (nbeginactual - nbegin*nbinsize)/(double) nbinsize;
	                             double dendfrac = ((nend+1)*nbinsize-nendactual-1)/(double) nbinsize;

				     if ((nbegin< posterior.length)&&(dbeginfrac>0))
				     {
				        //only removing if it would of made it into the actual posterior counts
				        float[] posterior_nbegin = posterior[nbegin];
			                for (int nstate = 0; nstate < numposteriorstates; nstate++)
			                {
		                           tallyoverlaplabel_ncoordfile[nstate]-=damount*dbeginfrac*posterior_nbegin[nstate];
					}
				     }			          
			          	
				     if ((nend < posterior.length)&&(dendfrac>0))
				     {
				         //only removing if it would of made it into the actual posterior counts		          
				        float[] posterior_nend = posterior[nend];
			                for (int nstate = 0; nstate < numposteriorstates; nstate++)
			                {
			                   tallyoverlaplabel_ncoordfile[nstate]-=damount*dendfrac*posterior_nend[nstate];
					}
				     }
				  }
			       }
			    }         
			 }
		      }		      
		      brcoords.close();		      
		   }
		}
	    }	 
	}

	if (!bposteriorfound)
	{
	    throw new IllegalArgumentException("no posterior file found for cell type "+szcell);
	}  

	for (int nfile = 0; nfile < dsumoverlaplabel.length; nfile++)
	{
	    double[] tallyoverlaplabel_nfile = tallyoverlaplabel[nfile];
	    for (int nindex = 0; nindex < tallyoverlaplabel_nfile.length; nindex++)
	    {
		dsumoverlaplabel[nfile] += tallyoverlaplabel_nfile[nindex];
	    }

	    if (dsumoverlaplabel[nfile] < 0.00000001)
	    {
	       throw new IllegalArgumentException("Coordinates in "+files[nfile]+" not assigned to any state");
	    }
	}



	outputenrichment(szoutfile, files,tallyoverlaplabel, tallylabel, dsumoverlaplabel,theColor,bcolscaleheat,ChromHMM.convertCharOrderToStringOrder(chorder),
			 sztitle,1,szlabelmapping,chorder); 

    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
     /**
      * Computes an enrichment based on a hard segmentation
      * szinputsegment - the file with the segmentation
      * szinputcoorddir - the directory containing the coordinates for enrichment
      * szinputcoordlist - if non-null only computes enrichments for file names in the specified
      * noffsetleft - the amount that should be subtracted from the left coordinate so it is 0-based inclusive
      * noffsetright - the amount that should be subtracted from the right coordinate so it is 0-based inclusive
      * nbinsize -  the binsize used in the segmentation
      * bcenter - if true only uses the center of the interval for determining enrichments
      * bunique - only counts a bin once, this is not applicable if bbaseres is true or busesignal is true
      * busesignal - if true then uses signal information about the interval if available
      * szcolfields - comma delimited list of the 0-based columns of chromosme,start,end, and if busesignal is true the signal column
      * bbaseres - if true gives enrichment at the base resolution opposed to the bin resolution
      * szoutfile - the file to which to output the 
      * bcolscaleheat - if true scales the heatmap individual for each column otherwise uses the same for all columns
      * theColor - theColor to use for the heatmap
      */
     public static void enrichmentMax(String szinputsegment,String szinputcoorddir,String szinputcoordlist,
				     int noffsetleft, int noffsetright,
                                     int nbinsize, boolean bcenter,boolean bunique, boolean busesignal,String szcolfields,
				      boolean bbaseres, String szoutfile,boolean bcolscaleheat,Color theColor,String sztitle, String szlabelmapping) throws IOException
    {
	
	ArrayList alsegments = new ArrayList(); //stores all the segments
	ArrayList alchromindex = new ArrayList();  //stores the index of the chromosome
	
	if (busesignal)
	{
	    bunique = false;
	}

	String szLine;
	HashMap hmchromMax = new HashMap(); //maps chromosome to the maximum index
	HashMap hmchromToIndex = new HashMap(); //maps chromosome to an index
	int nmaxlabel=0; // the maximum label found
	String szlabel="";
	//reads in the segmentation recording maximum position for each chromosome and
	//maximum label
	BufferedReader brinputsegment = Util.getBufferedReader(szinputsegment);
	while ((szLine = brinputsegment.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t ");
	    String szchrom = st.nextToken();
	    int nbegincoord = Integer.parseInt(st.nextToken());
	    int nendcoord = Integer.parseInt(st.nextToken());
	    if (nbegincoord % nbinsize != 0)
	    {
		throw new IllegalArgumentException("Binsize of "+nbinsize+" does not agree with input segment "+szLine);
	    }
	    int nbegin = nbegincoord/nbinsize;
	    int nend = (nendcoord-1)/nbinsize;
	    szlabel = st.nextToken();
	    short slabel;
            try
	    {
	       slabel  = (short) (Short.parseShort(szlabel));
	    }
	    catch (NumberFormatException ex)
	    {
               slabel  = (short) (Short.parseShort(szlabel.substring(1)));
	    }

	    alsegments.add(new SegmentRec(szchrom,nbegin,nend,slabel));
	    if (slabel > nmaxlabel)
	    {
		nmaxlabel = slabel;
	    }

	    Integer objMax = (Integer) hmchromMax.get(szchrom);
	    if (objMax == null)
	    {
		hmchromMax.put(szchrom,Integer.valueOf(nend));
		hmchromToIndex.put(szchrom, Integer.valueOf(hmchromToIndex.size()));
		alchromindex.add(szchrom);
	    }
	    else
	    {
		int ncurrmax = objMax.intValue();
		if (ncurrmax < nend)
		{
		    hmchromMax.put(szchrom, Integer.valueOf(nend));		    
		}
	    }
	}
	brinputsegment.close();

	int numchroms = alchromindex.size();
	short[][] labels = new short[numchroms][]; //stores the hard label assignments

	for (int nchrom = 0; nchrom < numchroms; nchrom++)
	{
	    int nsize = ((Integer) hmchromMax.get(alchromindex.get(nchrom))).intValue()+1;
	    labels[nchrom] = new short[nsize];
	    short[] labels_nchrom = labels[nchrom];
	    //sets to -1 so missing segments not counted as label 0
	    for (int npos = 0; npos < nsize; npos++)
	    {
		labels_nchrom[npos] = -1;
	    }
	       
	}	

	double[] tallylabel = new double[nmaxlabel+1];

	int numlines = alsegments.size();

	for (int nindex = 0; nindex < numlines; nindex++)
	{
	    SegmentRec theSegmentRec = (SegmentRec) alsegments.get(nindex);
	    int nbegin = theSegmentRec.nbegin;
	    int nend = theSegmentRec.nend;
	    short slabel = theSegmentRec.slabel;
	    int nchrom = ((Integer) hmchromToIndex.get(theSegmentRec.szchrom)).intValue();
	    short[] labels_nchrom = labels[nchrom];
	    //stores each label position in the genome
	    for (int npos = nbegin; npos <= nend; npos++)
	    {
		labels_nchrom[npos] = slabel;
		tallylabel[slabel]++; 
	    }
	}


	String[] files;

	if (szinputcoordlist == null)
        {
	    File dir = new File(szinputcoorddir);
	    //we don't have a specific list of files to include
	    //will use all files in the directory
	    if (dir.isDirectory())	  
	    {
		//throw new IllegalArgumentException(szinputcoorddir+" is not a directory!");
	   
	       files = dir.list(); 
	       Arrays.sort(files);
	       szinputcoorddir += "/";
	    }
	    else
	    {
		files = new String[1];
		files[0] = szinputcoorddir;
		szinputcoorddir = "";
	    }
	}
	else
	{
	    szinputcoorddir += "/";
	    //store in files all file names given in szinputcoordlist
	    BufferedReader brfiles = Util.getBufferedReader(szinputcoordlist);
	    ArrayList alfiles = new ArrayList();
	    while ((szLine = brfiles.readLine())!=null)
	    {
		alfiles.add(szLine);
	    }
	    brfiles.close(); 
	    files = new String[alfiles.size()];
	    for (int nfile = 0; nfile < files.length; nfile++)
	    {
		files[nfile] = (String) alfiles.get(nfile);
	    }
	}

	//for each enrichment category and state label gives a count of how often
	//overlapped by a segment optionally with signal
	double[][] tallyoverlaplabel = new double[files.length][nmaxlabel+1]; 
	double[] dsumoverlaplabel = new double[files.length];

	for (int nfile = 0; nfile < files.length; nfile++)
	{
	    double[] tallyoverlaplabel_nfile = tallyoverlaplabel[nfile];

	    int nchromindex = 0;
	    int nstartindex = 1;
	    int nendindex = 2;
	    int nsignalindex = 3;

	    if (szcolfields  != null)
	    {
	       StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
	       nchromindex = Integer.parseInt(stcolfields.nextToken());
	       nstartindex = Integer.parseInt(stcolfields.nextToken());
	       nendindex = Integer.parseInt(stcolfields.nextToken());

	       if (busesignal)
	       {
		   nsignalindex = Integer.parseInt(stcolfields.nextToken());
	       }
	    }


	    if (bunique)
	    {
	       Iterator itrChroms = hmchromToIndex.entrySet().iterator();
	       while (itrChroms.hasNext())
	       {
	          Map.Entry pairs = (Map.Entry) itrChroms.next();
		  String szchrom =(String) pairs.getKey();
		  int nchrom = ((Integer) pairs.getValue()).intValue();
	          short[] labels_nchrom = labels[nchrom];

	          //reading in the coordinates to overlap with
                  BufferedReader brcoords = Util.getBufferedReader(szinputcoorddir +files[nfile]);
		  ArrayList alrecs = new ArrayList();
	          while ((szLine = brcoords.readLine())!=null)
	          {
	             if (szLine.trim().equals("")) continue;
	             String[] szLineA = szLine.split("\\s+");
		     if (nstartindex >= szLineA.length)
		     {
			 throw new IllegalArgumentException(nstartindex+" is an invalid index for "+szLine+" in "+szinputcoorddir+files[nfile]);
		     }


                     if (nendindex >= szLineA.length)
		     {
			 throw new IllegalArgumentException(nendindex+" is an invalid index for "+szLine+" in "+szinputcoorddir+files[nfile]);
		     }


	             String szcurrchrom = szLineA[nchromindex];
		     if (szchrom.equals(szcurrchrom))
		     {
		        int nbeginactual =Integer.parseInt(szLineA[nstartindex])-noffsetleft;
		        int nendactual =Integer.parseInt(szLineA[nendindex])-noffsetright;
			if (bcenter)
			{
			    nbeginactual = (nbeginactual+nendactual)/2;
			    nendactual = nbeginactual;
			}
			alrecs.add(new Interval(nbeginactual,nendactual));
		     }
		  }
		  brcoords.close();

		  Object[] alrecA = alrecs.toArray();
		  Arrays.sort(alrecA,new IntervalCompare());

		  boolean bclosed = true;
		  int nintervalstart = -1;
		  int nintervalend = -1;
		  boolean bdone = false;

		  for (int nindex = 0; (nindex <= alrecA.length&&(alrecA.length>0)); nindex++)
		  {
		      int ncurrstart=-1;
		      int ncurrend=-1;


		     if (nindex == alrecA.length)
		     {
			 bdone = true;
		     }
		     else 
		     {
		        ncurrstart = ((Interval) alrecA[nindex]).nstart;
		        ncurrend = ((Interval) alrecA[nindex]).nend;
                        if (nindex == 0)
		        {
			   nintervalstart = ncurrstart;
			   nintervalend = ncurrend;
		        }
		        else if (ncurrstart <= nintervalend)
		        {
			    //this read is still in the active interval
			    //extending the current active interval 
		           if (ncurrend > nintervalend)
		           {
		              nintervalend = ncurrend;
		           }
			}		     
		        else 
		        {
		           //just finished the current active interval
			   bdone = true;
			}
		     }

		     if (bdone)
		     {		      						
		        int nbegin = nintervalstart/nbinsize;
		        int nend = nintervalend/nbinsize;
			if (nend >= labels_nchrom.length)
			{
			   nend = labels_nchrom.length - 1;
			}

	                for (int nbin = nbegin; nbin <= nend; nbin++)
	                {
			   if (labels_nchrom[nbin]>=0)
			   {
			       tallyoverlaplabel_nfile[labels_nchrom[nbin]]++;
			   }
			}

		        if (bbaseres)
		        { 
		           //dbeginfrac represents the fraction of bases the nbegin interval
		           //which came after the actual nbeginactual
	                   double dbeginfrac = (nintervalstart - nbegin*nbinsize)/(double) nbinsize;
			   
			   //dendfrac represents the fraction of bases after the end position in the interval
	                   double dendfrac = ((nend+1)*nbinsize-nintervalend-1)/(double) nbinsize;
			   
			   if ((nbegin < labels_nchrom.length)&&(labels_nchrom[nbegin]>=0)&&(dbeginfrac>0))
		           { 
			      
		              //only counted the bases if nbegin was less than labels_nchrom.length  
		              tallyoverlaplabel_nfile[labels_nchrom[nbegin]]-=dbeginfrac;
		           }

                           if ((nend < labels_nchrom.length)&&(labels_nchrom[nend]>=0)&&(dendfrac>0))
		           {
		              //only counted the bases if nend was less than labels_nchrom.length  
		              tallyoverlaplabel_nfile[labels_nchrom[nend]]-=dendfrac;
		           }
			}			   

		        nintervalstart = ncurrstart; 
		        nintervalend = ncurrend;
			bdone = false;
		     }		  
		  }
	       }
	    }
	    else
	    {
	       BufferedReader brcoords = Util.getBufferedReader(szinputcoorddir +files[nfile]);
	       while ((szLine = brcoords.readLine())!=null)
	       {
	          if (szLine.trim().equals("")) continue;
	          String[] szLineA = szLine.split("\\s+");

	          String szchrom = szLineA[nchromindex];
	          int nbeginactual =Integer.parseInt(szLineA[nstartindex])-noffsetleft;
	          int nbegin = nbeginactual/nbinsize;

		  int nendactual =Integer.parseInt(szLineA[nendindex])-noffsetright;
		  int nend = nendactual/nbinsize;

		  double damount;
	          if ((busesignal)&&(nsignalindex < szLineA.length))
	          {
	      	     damount = Double.parseDouble(szLineA[nsignalindex]);
		  }
	          else
	          {
	             damount = 1;
	          }

	          Integer objChrom = (Integer) hmchromToIndex.get(szchrom);
	          if (objChrom != null)
	          {
		     //we have the chromosome corresponding to this read
	             int nchrom = objChrom.intValue();
		     short[] labels_nchrom = labels[nchrom];

		     if (bcenter)
	             {
		         //using the center position of the interval only
		        int ncenter = (nbeginactual+nendactual)/(2*nbinsize);
		        if ((ncenter < labels_nchrom.length)&&(labels_nchrom[ncenter]>=0))
		        {
		           tallyoverlaplabel_nfile[labels_nchrom[ncenter]]+=damount;			   
		        }
		     }
	             else
	             {
		        //using the full interval range
			//no requirement on uniqueness
		        if (nend >= labels_nchrom.length)
			{
			   nend = labels_nchrom.length - 1;
			}

			
	                for (int nindex = nbegin; nindex <= nend; nindex++)
	                {
			   if (labels_nchrom[nindex]>=0)
			   {
			      //increment overlap tally not checking for uniqueness
 	                      tallyoverlaplabel_nfile[labels_nchrom[nindex]]+=damount;
			   }
			}	       

			if (bbaseres)
			{ 
			   //dbeginfrac represents the fraction of bases the nbegin interval
			   //which came after the actual nbeginactual
	                   double dbeginfrac = (nbeginactual - nbegin*nbinsize)/(double) nbinsize;
			   
			   //dendfrac represents the fraction of bases after the end position in the interval
	                   double dendfrac = ((nend+1)*nbinsize-nendactual-1)/(double) nbinsize;

			   if ((nbegin < labels_nchrom.length)&&(labels_nchrom[nbegin]>=0)&&(dbeginfrac>0))
			   { 
			      //only counted the bases if nbegin was less than labels_nchrom.length  
			      tallyoverlaplabel_nfile[labels_nchrom[nbegin]]-=damount*dbeginfrac;
			   }

                           if ((nend < labels_nchrom.length)&&(labels_nchrom[nend]>=0)&&(dendfrac>0))
		           {
			      //only counted the bases if nend was less than labels_nchrom.length  
			      tallyoverlaplabel_nfile[labels_nchrom[nend]]-=damount*dendfrac;
			   }			   
			}
		     }
		  }
	       }
	       brcoords.close();
	    }

	    for (int nindex = 0; nindex < tallyoverlaplabel_nfile.length; nindex++)
	    {
		dsumoverlaplabel[nfile] += tallyoverlaplabel_nfile[nindex];
	    }
		
            if (dsumoverlaplabel[nfile] < 0.00000001)
	    {
	       throw new IllegalArgumentException("Coordinates in "+files[nfile]+" not assigned to any state");
	    }
	}

	outputenrichment(szoutfile, files,tallyoverlaplabel, tallylabel, dsumoverlaplabel,theColor,
			 bcolscaleheat,ChromHMM.convertCharOrderToStringOrder(szlabel.charAt(0)),sztitle,0,szlabelmapping,szlabel.charAt(0));
    }

    /////////////////////////////////////////////////////////////////
    /**
     * Loads contents of szlabelmapping into a HashMap
     */
    private static HashMap makeLabelMapping(String szlabelmapping) throws IOException
    {
	HashMap hmlabelExtend = new HashMap();
	if (szlabelmapping != null)
	{
           BufferedReader bridlabel =  Util.getBufferedReader(szlabelmapping);
	   String szLine;

	   //Loading in a mapping from state ID to a label description

	   while ((szLine = bridlabel.readLine())!=null)
           {
	      StringTokenizer st = new StringTokenizer(szLine,"\t");
	      String szID = st.nextToken();
	      String szLabelExtend = st.nextToken();
	      hmlabelExtend.put(szID,szLabelExtend);
	   }
	   bridlabel.close();
	}

	return hmlabelExtend;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Prints the enrichment to a text file and heatmap for the enrichment tallies
     */
    private static void outputenrichment(String szoutfile,String[] files,double[][] tallyoverlaplabel, 
                                         double[] tallylabel, double[] dsumoverlaplabel,Color theColor,
                                         boolean bcolscaleheat,String szstateorder,String sztitle, int noffset,String szlabelmapping,char chorder) throws IOException
     {

	HashMap hmlabelExtend = makeLabelMapping(szlabelmapping);

	double dsumlabel = 0;
	for (int ni = 0; ni < tallylabel.length; ni++)
        {
	    dsumlabel += tallylabel[ni];
	}

        NumberFormat nf5 = NumberFormat.getInstance();
        nf5.setMaximumFractionDigits(5);
	nf5.setGroupingUsed(false);
	nf5.setMinimumFractionDigits(5);

	System.out.println("Writing to file "+szoutfile+".txt");
	PrintWriter pw = new PrintWriter(new FileWriter(szoutfile+".txt"));

	pw.print("state ("+szstateorder+" order)\tGenome %");
	for (int nfile = 0; nfile < files.length; nfile++)
	{
	    pw.print("\t"+files[nfile]);
	}
	pw.println();


	double[][] heatmapfold;


	if (bcolscaleheat)
	{
	    heatmapfold = new double[tallyoverlaplabel[0].length][tallyoverlaplabel.length+1];
	}
	else
	{
	    heatmapfold = new double[tallyoverlaplabel[0].length][tallyoverlaplabel.length];
	}

	int numelim = 0;

	for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	{
	    if (tallylabel[nstate] > 0)
	    {
		pw.print(nstate+noffset);
                String szsuffix;

                if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(nstate+noffset)))!=null)
		{
	      	   pw.print("_"+szsuffix);
	        }

		if (bcolscaleheat)
		{
		    //only include genome % if scaling by column
		   heatmapfold[nstate][0] =100*(tallylabel[nstate]/dsumlabel); //first column is heatmap
		   pw.print("\t"+nf5.format(heatmapfold[nstate][0]));
		}
		for (int nfile = 0; nfile < tallyoverlaplabel.length; nfile++)
	        {
		    //numerator is fraction of total of category in state while denominator is total fraction state is
	           double dfold = (tallyoverlaplabel[nfile][nstate]/dsumoverlaplabel[nfile])/
	                           (tallylabel[nstate]/dsumlabel);
		   if (bcolscaleheat)
		   {
		      heatmapfold[nstate][nfile+1] = dfold;
		   }
		   else
		   {
		      heatmapfold[nstate][nfile] = dfold;
		   }
		   pw.print("\t"+nf5.format(dfold));
		}
		pw.println();
	    }
	    else
	    {
		numelim++;
	    }
	}
	pw.print("Base\t100");
        for (int nfile = 0; nfile < tallyoverlaplabel.length; nfile++)
	{
	    pw.print("\t"+nf5.format(100*(dsumoverlaplabel[nfile]/dsumlabel)));
	}
	pw.println();
	pw.close();

	String[] rowlabels = new String[tallyoverlaplabel[0].length];
	if (numelim > 0)
	{
	   double[][] heatmapreduce = new double[heatmapfold.length-numelim][heatmapfold[0].length];
	   rowlabels = new String[heatmapreduce.length];
	   int nkeepindex = 0;
	   for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	   {
              if (tallylabel[nstate] > 0)
	      {
		  for (int ncol = 0; ncol < heatmapfold[nstate].length; ncol++)
		  {
		      heatmapreduce[nkeepindex][ncol] = heatmapfold[nstate][ncol];
		  }
		  rowlabels[nkeepindex] = ""+(nstate+noffset);
		  nkeepindex++;
	      }
	   }
	   heatmapfold = heatmapreduce;
	}
	else
	{
	    rowlabels = new String[tallyoverlaplabel[0].length];
	    for (int ni = 0; ni < rowlabels.length; ni++)
	    {
		rowlabels[ni] = ""+(ni+noffset);
	    }
	}

	for (int ni = 0; ni < rowlabels.length; ni++)
	{
	   String szsuffix;
	   if ((szsuffix = (String) hmlabelExtend.get(""+chorder+rowlabels[ni]))!=null)
           {
       	      rowlabels[ni] = rowlabels[ni]+"_"+szsuffix;
	   }
	}

	if (bcolscaleheat)
	{
	   for (int nfile = 0; nfile < heatmapfold[0].length; nfile++)
	   {	    
	      double dmaxval = Double.NEGATIVE_INFINITY;
	      double dminval = Double.POSITIVE_INFINITY;

	      for (int nstate = 0; nstate < heatmapfold.length; nstate++)
	      {
	         if (heatmapfold[nstate][nfile] > dmaxval)
		 {
		    dmaxval = heatmapfold[nstate][nfile];
		 }

	         if (heatmapfold[nstate][nfile] < dminval)
		 {
		    dminval = heatmapfold[nstate][nfile];
		 }
	      }
	
	      for (int nstate = 0; nstate < heatmapfold.length; nstate++)
	      {
		  heatmapfold[nstate][nfile] =(heatmapfold[nstate][nfile]-dminval)/dmaxval;
	      }
	   }
	}

	String[] collabels = new String[heatmapfold[0].length];
	if (bcolscaleheat)
	{
	   collabels[0] = "Genome %";
	}

	for (int ni = 0; ni < files.length; ni++)
	{
	    int nfirstindex = files[ni].indexOf('.');
	    if (nfirstindex == -1)
	    {
	       if (bcolscaleheat)
	       {
	          collabels[ni+1] = files[ni];
	       }
	       else
	       {
	          collabels[ni] = files[ni];
	       }
	    }
	    else
	    {
	       String szsub =  files[ni].substring(0,nfirstindex);

	       if (bcolscaleheat)
	       {
	          collabels[ni+1] = szsub;
	       }
	       else
	       {
	          collabels[ni] = szsub;
	       }
	    }
	}


	makeEnrichmentHeatMap(heatmapfold, collabels, rowlabels,szoutfile,theColor,sztitle,
			      "Category","State ("+szstateorder+" order)");
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Record have the index of the chromsome, coordinate position, strand, and signal
     */
    private static class RecAnchorIndex
    {
	int nchromindex;
        int npositionindex;
        int nstrandindex;
        int nsignalindex;

	RecAnchorIndex(int nchromindex, int npositionindex, int nstrandindex, int nsignalindex)
	{
	    this.nchromindex = nchromindex;
	    this.npositionindex = npositionindex;
	    this.nstrandindex = nstrandindex;
	    this.nsignalindex = nsignalindex;
	}
    }

    /**
     * Returns a Record of the indicies as determined by the contents of szcolfields, busestrand, and busesignal
     * of the chromosome, coordinate, strand and signal positions
     */
    private static RecAnchorIndex getAnchorIndex(String szcolfields, boolean busestrand, boolean busesignal)
    {
	int nchromindex=0;
        int npositionindex=1;
        int nstrandindex=2;
        int nsignalindex=3;

	if (busestrand && busesignal)
        {
	   if (szcolfields != null)
           {
	       StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
	       nchromindex = Integer.parseInt(stcolfields.nextToken());
	       npositionindex = Integer.parseInt(stcolfields.nextToken());
	       nstrandindex = Integer.parseInt(stcolfields.nextToken());
	       nsignalindex = Integer.parseInt(stcolfields.nextToken());
	   }
	}
        else if (busestrand && !busesignal)
        {
	   if (szcolfields != null)
           {
	       StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
	       nchromindex = Integer.parseInt(stcolfields.nextToken());
	       npositionindex = Integer.parseInt(stcolfields.nextToken());
	       nstrandindex = Integer.parseInt(stcolfields.nextToken());
	   }
	}
        else if (!busestrand && busesignal)
        {
	   if (szcolfields != null)
	   {
	      StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
	      nchromindex = Integer.parseInt(stcolfields.nextToken());
	      npositionindex = Integer.parseInt(stcolfields.nextToken());
	      nsignalindex = Integer.parseInt(stcolfields.nextToken());
	   }
	   else if (szcolfields == null)
           {
	       nsignalindex = 2;
	   }
	}
        else //both are off
	{
	   if (szcolfields != null)
	   {
	      StringTokenizer stcolfields = new StringTokenizer(szcolfields,",");
	      nchromindex = Integer.parseInt(stcolfields.nextToken());
	      npositionindex = Integer.parseInt(stcolfields.nextToken());
	   }
	}

	return new RecAnchorIndex(nchromindex, npositionindex, nstrandindex, nsignalindex);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Outputs a neighborhood enrichment based on a given hard segmentation
     * szinputsegmentation - is the input segmentation ; assumes segments are non-negative integers
     * szanchorpositions - is a file with the coordinates of the anchor positions
     * nbinsize - the binsize used in the segmentation, this ensures any position at this bin size 
     * numleft -  the number of intervals to the left of the anchor position to display the heatmap for
     * numright - the number of intervals to the right of the anchor position to display the heatmap for
     * nspacing - the frequency at which to display enrichments
     * noffsetanchor - the amount that should be substract so the anchor coordinates are 0 based
     * busestrand - true if strand orientation associated with the position is given
     * busesignal - true if signal associated with the data file should be used
     * szcolfields - a comma list of the 0-based index of 
     * if (busestrand and busesignal)     chromosome,pos,strand,signal default is 0,1,2,3
     * if (busestrand and not busesignal) chromosome,pos,strand default is 0,1,2
     * if (not busestrand and busesignal) chromosome,pos,signal default is 0,1,2
     * if (not busestrand not busesignal) chromosome,pos default is 0,1
     * szoutfile - the name of the text file for which the enrichments should be displayed
     * heatmaps are written to the same place with the extensions '.png' and '.svg'
     * theColor - theColor to use for the heatmap program
     */ 
     public static void neighborhoodMax(String szinputsegmentation,String szanchorpositions,
                                       int nbinsize, int numleft, int numright, int nspacing, 
					boolean busestrand, boolean busesignal, String szcolfields,
					int noffsetanchor, String szoutfile,Color theColor, 
					String sztitle,String szlabelmapping) throws IOException
    {
	//stores all the segments in the data
	ArrayList alsegments = new ArrayList();

	//an array of chromosome names
	ArrayList alchromindex = new ArrayList();

	String szLine;

	//stores the largest index value for each chromosome
	HashMap hmchromMax = new HashMap();

	//maps chromosome names to index values
	HashMap hmchromToIndex = new HashMap();

	//stores the maximum integer label value
	int nmaxlabel=0;
	String szlabel ="";
	BufferedReader brinputsegment = Util.getBufferedReader(szinputsegmentation);

	//this loops reads in the segmentation 
	while ((szLine = brinputsegment.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t ");
	    String szchrom = st.nextToken();
            //assumes segments are in standard bed format which to get to 
	    //0-based inclusive requires substract 1 from the end
	    int nbegin = Integer.parseInt(st.nextToken())/nbinsize;
	    int nend = (Integer.parseInt(st.nextToken())-1)/nbinsize; 
	    szlabel = st.nextToken();
	    short slabel;
            try
	    {
		slabel  = (short) (Short.parseShort(szlabel));
	    }
	    catch (NumberFormatException ex)
	    {
               slabel  = (short) (Short.parseShort(szlabel.substring(1)));
	    }

	    alsegments.add(new SegmentRec(szchrom,nbegin,nend,slabel));

	    if (slabel > nmaxlabel)
	    {
		nmaxlabel = slabel;
	    }

	    Integer objMax = (Integer) hmchromMax.get(szchrom);
	    if (objMax == null)
	    {
		hmchromMax.put(szchrom,Integer.valueOf(nend));
		hmchromToIndex.put(szchrom, Integer.valueOf(hmchromToIndex.size()));
		alchromindex.add(szchrom);
	    }
	    else
	    {
		int ncurrmax = objMax.intValue();
		if (ncurrmax < nend)
		{
		    hmchromMax.put(szchrom, Integer.valueOf(nend));		    
		}
	    }
	}
	brinputsegment.close();

	int numchroms = alchromindex.size();
        short[][] labels = new short[numchroms][];
	//allocates space store all the segment labels for each chromosome
	for (int nchrom = 0; nchrom < numchroms; nchrom++)
	{
	    int nsize = ((Integer) hmchromMax.get(alchromindex.get(nchrom))).intValue()+1;
	    labels[nchrom] = new short[nsize];
	    short[] labels_nchrom = labels[nchrom];
	    for (int npos = 0; npos < nsize; npos++)
	    {
		labels_nchrom[npos] = -1;
	    }
	}	

	//a tally on how frequently each label occurs
	double[] tallylabel = new double[nmaxlabel+1];

	int numlines = alsegments.size();

	//this loop stores into labels the full segmentation
	//and a count of how often each label occurs
	for (int nindex = 0; nindex < numlines; nindex++)
	{
	    SegmentRec theSegmentRec = (SegmentRec) alsegments.get(nindex);
	    int nchrom = ((Integer) hmchromToIndex.get(theSegmentRec.szchrom)).intValue();
	    short[] labels_nchrom = labels[nchrom];
	    int nbegin = theSegmentRec.nbegin;
	    int nend  = theSegmentRec.nend;
	    short slabel = theSegmentRec.slabel;
	    for (int npos = nbegin; npos <= nend; npos++)
	    {
		labels_nchrom[npos] = slabel;
		if (slabel >= 0)
		{
		   tallylabel[slabel]++; 
		}
	    }
	}

        //the number of additional intervals to the left and right to include

	//the center anchor position
	int numintervals = 1+numleft+numright;

	//stores a tally for each position relative to an anchor how frequently the label was observed
	double[][] tallyoverlaplabel = new double[numintervals][nmaxlabel+1]; 

	//stores a tally for the total signal associated with each anchor position
	double[] dsumoverlaplabel = new double[numintervals];

	RecAnchorIndex theAnchorIndex = getAnchorIndex(szcolfields, busestrand, busesignal);

	//reads in the anchor position 
        BufferedReader brcoords = Util.getBufferedReader(szanchorpositions);
	while ((szLine = brcoords.readLine())!=null)
        {
	    if (szLine.trim().equals("")) continue;
	   String[] szLineA = szLine.split("\\s+");

           String szchrom = szLineA[theAnchorIndex.nchromindex];        
	   int nanchor = (Integer.parseInt(szLineA[theAnchorIndex.npositionindex])-noffsetanchor);
	   boolean bposstrand = true;
	   if (busestrand)
	   {
	      String szstrand = szLineA[theAnchorIndex.nstrandindex];	    
	      if (szstrand.equals("+"))
	      {
	         bposstrand = true;
	      }
              else if (szstrand.equals("-"))
              {
      	         bposstrand = false;
	      }
	      else
	      {
      	         throw new IllegalArgumentException(szstrand +" is an invalid strand. Strand should be '+' or '-'");
	      }	      
	   }

	   double damount;

	   if ((busesignal)&&(theAnchorIndex.nsignalindex< szLineA.length))
	   {
	       damount = Double.parseDouble(szLineA[theAnchorIndex.nsignalindex]);
	   }
	   else
           {
	      damount = 1;
	   }

	   //updates the tallys for the given anchor position
	   Integer objChrom = (Integer) hmchromToIndex.get(szchrom);
           if (objChrom != null)
	   {
	      int nchrom = objChrom.intValue();
	      short[] labels_nchrom = labels[nchrom];

	      if (bposstrand)
	      {
	         int ntallyindex = 0;
	         for(int noffset= -numleft; noffset <= numright; noffset++)
	         {
		    int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		    if ((nposindex >=0)&&(nposindex < labels_nchrom.length)&&(labels_nchrom[nposindex]>=0))
		    {
	               tallyoverlaplabel[ntallyindex][labels_nchrom[nposindex]] += damount;		      
		    }
		    ntallyindex++;
		 }
	      }
	      else
	      {
	         int ntallyindex = 0;
	         for(int noffset= numright; noffset >= -numleft; noffset--)
	         {
		     int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		     if ((nposindex >=0)&&(nposindex < labels_nchrom.length)&&(labels_nchrom[nposindex]>=0))
		    {
	               tallyoverlaplabel[ntallyindex][labels_nchrom[nposindex]]+=damount;		      
		    }
		    ntallyindex++;
		 }
	      }
	   }
	}
        brcoords.close(); 	    

	outputneighborhood(tallyoverlaplabel,tallylabel,dsumoverlaplabel,szoutfile,nspacing,numright,
                           numleft,theColor,ChromHMM.convertCharOrderToStringOrder(szlabel.charAt(0)),sztitle,0,szlabelmapping,szlabel.charAt(0));
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Outputs a signal enrichment agreggation 
     * szinputsegmentation - szposteriordir the directory with the posterior files
     * szcell - the cell type to use for the neighborhood enrichments
     * szanchorpositions - is a file with the coordinates of the anchor positions
     * nbinsize - the binsize used in the segmentation, this ensures any position at this bin size 
     * numleft -  the number of intervals to the left of the anchor position to display the heatmap for
     * numright - the number of intervals to the right of the anchor position to display the heatmap for
     * nspacing - the frequency at which to display enrichments
     * noffsetanchor - the amount that should be substract so the anchor coordinates are 0 based
     * busestrand - true if strand orientation associated with the position is given
     * busesignal - true if signal associated with the data file should be used
     * szcolfields - a comma list of the 0-based index of 
     * if (busestrand and busesignal)     chromosome,pos,strand,signal default is 0,1,2,3
     * if (busestrand and not busesignal) chromosome,pos,strand default is 0,1,2
     * if (not busestrand and busesignal) chromosome,pos,signal default is 0,1,2
     * if (not busestrand not busesignal) chromosome,pos default is 0,1
     * szoutfile - the name of the text file for which the enrichments should be displayed
     * a heatmap is written to the same place with the extensions '.png' and '.svg'
     * theColor - theColor to use for the heatmap program
     */
     public static void neighborhoodSignal(String szposteriordir,String szcell,String szanchorpositions, int nbinsize, int numleft, int numright, int nspacing,
					     boolean busestrand, boolean busesignal,String szcolfields, int noffsetanchor,
					   String szoutfile,Color theColor,String sztitle,String szlabelmapping) throws IOException
    {
	//posterior here is really signal just using equivalent variable names
	//list of possible posterior files
	File posteriordir = new File(szposteriordir);
	String[] posteriorfiles = posteriordir.list();

	String szLine;
	int numlocs = 0;

	int numposteriorstates = 0;
	char chorder= '0';//just need to initalize it to something
	int nfirstindex = 0; //the index of the first file found
	boolean bfirst = true;

	int numintervals = 1+numleft+numright; //total number of intervals

	double[][] tallyoverlaplabel = null; //tally of the overlap of each interval and state
	double[] tallylabel = null; //tally of the overlap of each interval

	double[] dsumoverlaplabel = new double[numintervals];
	String szmarknames="";
	boolean bposteriorfound = false;
	for (int nfile = 0; nfile < posteriorfiles.length; nfile++)
	{	
	    String szposteriorfiles_nfile = posteriorfiles[nfile];

	    //going through all posterior files
	    if (szposteriorfiles_nfile.contains("_signal"))
	    {		
		BufferedReader brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		szLine = brposterior.readLine();
	        if (szLine ==null)
	        {
	           throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" is empty!");
	        }
		StringTokenizer st =new StringTokenizer(szLine,"\t");
		String szcurrcell = st.nextToken();
		if ((!szcurrcell.equals(szcell))&&(!szcell.equals("")))
		{
		    brposterior.close();
		}
		else
		{
		   bposteriorfound = true;

  	           //must match cell type or consistent with empty cell type

		   String szchrom = st.nextToken();
	    
		   int numlines = 0;
		   szLine = brposterior.readLine(); //gets state header
		   if (szLine == null)
		   {
		       throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" only has one line!");
		   }
		   st = new StringTokenizer(szLine,"\t");
		   int numcurrstates = st.countTokens();
	           if (bfirst)
	           {
		       chorder = st.nextToken().charAt(0);
		       bfirst = false;
		       nfirstindex = nfile;
		       numposteriorstates = numcurrstates;
		       tallylabel = new double[numposteriorstates];
		       tallyoverlaplabel = new double[numintervals][numposteriorstates];
		   }
	           else if (numposteriorstates != numcurrstates)
	           {
	              throw new IllegalArgumentException("Number of states "+numcurrstates+" in "+szposteriorfiles_nfile+
							" does not match number of states "+numposteriorstates+" in "+posteriorfiles[nfirstindex]);
		   }

	           while ((szLine = brposterior.readLine())!=null)
	           {
		      numlines++;
		   }
		   brposterior.close();	     
		   numlocs += numlines;
	           //loading in the posterior data
	           float[][] posterior = new float[numlines][numposteriorstates];
		
		   brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		   int nline = 0;
	           brposterior.readLine();//specifies the header
	           szmarknames = brposterior.readLine();
	           while ((szLine = brposterior.readLine())!=null)
	           {
		      st = new StringTokenizer(szLine,"\t ");
		      float[] posterior_nline = posterior[nline];
                      for (int nstate = 0; nstate < numposteriorstates; nstate++)
                      {
		         float fval = Float.parseFloat(st.nextToken());
		         posterior_nline[nstate] = fval;
		         //keeps a tally of often each state occurs probabilistically
		         tallylabel[nstate] += fval;
		      }
		      nline++;
		   }
		   brposterior.close();

		   boolean[] counted = null;
 
		   BufferedReader brcoords = Util.getBufferedReader(szanchorpositions);
		   RecAnchorIndex theAnchorIndex = getAnchorIndex(szcolfields, busestrand, busesignal);
	           while ((szLine = brcoords.readLine())!=null)
	           {
		      if (szLine.trim().equals("")) continue;
	              //going through all the coordinates keeping only those for the current chromosome
	              String[] szLineA = szLine.split("\\s+");
                      String szreadchrom = szLineA[theAnchorIndex.nchromindex];        

	              if (szreadchrom.equals(szchrom))
	              {  
	                 int nanchor = (Integer.parseInt(szLineA[theAnchorIndex.npositionindex])-noffsetanchor);
	                 boolean bposstrand = true;
	                 if (busestrand)
	                 {
	                    String szstrand = szLineA[theAnchorIndex.nstrandindex];	    
	                    if (szstrand.equals("+"))
	                    {
			       bposstrand = true;
			    }
                            else if (szstrand.equals("-"))
                            {
			       bposstrand = false;
			    }
	                    else
	                    {
			       throw new IllegalArgumentException(szstrand +" is an valid strand. Strand should be '+' or '-'");
			    }     
			 }
		 
		         double damount;

	                 if ((busesignal)&&(theAnchorIndex.nsignalindex < szLineA.length))
	                 {
	                    damount = Double.parseDouble(szLineA[theAnchorIndex.nsignalindex]);
			 }
	                 else
                         {
	                    damount = 1;
			 }

	                 if (bposstrand)
	                 {
		            int ntallyindex = 0;
	                    for(int noffset= -numleft; noffset <= numright; noffset++)
	                    {
		               int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		               if ((nposindex >=0)&&(nposindex < posterior.length))
		               {
			          //within interval
		                  float[] posterior_nposindex = posterior[nposindex];
		                  double[] tallyoverlaplabel_ntallyindex = tallyoverlaplabel[ntallyindex];
		                  for (int nstate = 0; nstate < numposteriorstates; nstate++)
		                  {
				     //increments each state by the posterior amount
		                     tallyoverlaplabel_ntallyindex[nstate]+= damount*posterior_nposindex[nstate];		      
				  }
				  dsumoverlaplabel[ntallyindex]++;
			       }
			       ntallyindex++;
			    }
			 }
	                 else
	                 {
		            int ntallyindex = 0;
	                    for(int noffset= numright; noffset >= -numleft; noffset--)
	                    {
		               int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		               if ((nposindex >=0)&&(nposindex < posterior.length))
		               {
			          float[] posterior_nposindex = posterior[nposindex];
			          double[] tallyoverlaplabel_ntallyindex = tallyoverlaplabel[ntallyindex];
		                  for (int nstate = 0; nstate < numposteriorstates; nstate++)
		                  {
				     //increments each state by the posterior amount
	                             tallyoverlaplabel_ntallyindex[nstate]+=damount*posterior_nposindex[nstate];		      
				  }
				  dsumoverlaplabel[ntallyindex]++;
			       }
			       ntallyindex++;		    			   
			    }
			 }
		      }
		   }		
		   brcoords.close();
		}
	    }
	}	

	if (!bposteriorfound)
	{
	    throw new IllegalArgumentException("No posterior file found for cell type "+szcell);
	}   
 
	outputneighborhoodsignal(tallyoverlaplabel,tallylabel,dsumoverlaplabel,numlocs,szoutfile,nspacing,numright,
				 numleft,theColor,ChromHMM.convertCharOrderToStringOrder(chorder),sztitle,szmarknames,szlabelmapping,chorder);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Outputs a neighborhood enrichment based on a given posterior segmentaation
     * szinputsegmentation - szposteriordir the directory with the posterior files
     * szcell - the cell type to use for the neighborhood enrichments
     * szanchorpositions - is a file with the coordinates of the anchor positions
     * nbinsize - the binsize used in the segmentation, this ensures any position at this bin size 
     * numleft -  the number of intervals to the left of the anchor position to display the heatmap for
     * numright - the number of intervals to the right of the anchor position to display the heatmap for
     * nspacing - the frequency at which to display enrichments
     * noffsetanchor - the amount that should be substract so the anchor coordinates are 0 based
     * busestrand - true if strand orientation associated with the position is given
     * busesignal - true if signal associated with the data file should be used
     * szcolfields - a comma list of the 0-based index of 
     * if (busestrand and busesignal)     chromosome,pos,strand,signal default is 0,1,2,3
     * if (busestrand and not busesignal) chromosome,pos,strand default is 0,1,2
     * if (not busestrand and busesignal) chromosome,pos,signal default is 0,1,2
     * if (not busestrand not busesignal) chromosome,pos default is 0,1
     * szoutfile - the name of the text file for which the enrichments should be displayed
     * a heatmap is written to the same place with the extensions '.png' and '.svg'
     * theColor - theColor to use for the heatmap program
     */
     public static void neighborhoodPosterior(String szposteriordir,String szcell,String szanchorpositions, int nbinsize, int numleft, int numright, int nspacing,
					     boolean busestrand, boolean busesignal,String szcolfields, int noffsetanchor,
					      String szoutfile,Color theColor,String sztitle,String szlabelmapping) throws IOException
    {

	//list of possible posterior files
	File posteriordir = new File(szposteriordir);
	String[] posteriorfiles = posteriordir.list();

	String szLine;

	int numposteriorstates = 0;
	char chorder= '0';//just need to initalize it to something
	int nfirstindex = 0; //the index of the first file found
	boolean bfirst = true;

	int numintervals = 1+numleft+numright; //total number of intervals

	double[][] tallyoverlaplabel = null; //tally of the overlap of each interval and state
	double[] tallylabel = null; //tally of the overlap of each interval

	double[] dsumoverlaplabel = new double[numintervals];

	boolean bposteriorfound = false;
	for (int nfile = 0; nfile < posteriorfiles.length; nfile++)
	{	
	    String szposteriorfiles_nfile = posteriorfiles[nfile];

	    //going through all posterior files
	    if (szposteriorfiles_nfile.contains("_posterior"))
	    {		
		BufferedReader brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		szLine = brposterior.readLine();
	        if (szLine ==null)
	        {
	           throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" is empty!");
	        }
		StringTokenizer st =new StringTokenizer(szLine,"\t");
		String szcurrcell = st.nextToken();
		if ((!szcurrcell.equals(szcell))&&(!szcell.equals("")))
		{
		    brposterior.close();
		}
		else
		{
		   bposteriorfound = true;

  	           //must match cell type or consistent with empty cell type

		   String szchrom = st.nextToken();
	    
		   int numlines = 0;
		   szLine = brposterior.readLine(); //gets state header
		   if (szLine == null)
		   {
		       throw new IllegalArgumentException(szposteriordir+"/"+szposteriorfiles_nfile+" only has one line!");
		   }
		   st = new StringTokenizer(szLine,"\t");
		   int numcurrstates = st.countTokens();
	           if (bfirst)
	           {
		       chorder = st.nextToken().charAt(0);
		       bfirst = false;
		       nfirstindex = nfile;
		       numposteriorstates = numcurrstates;
		       tallylabel = new double[numposteriorstates];
		       tallyoverlaplabel = new double[numintervals][numposteriorstates];
		   }
	           else if (numposteriorstates != numcurrstates)
	           {
	              throw new IllegalArgumentException("Number of states "+numcurrstates+" in "+szposteriorfiles_nfile+
							" does not match number of states "+numposteriorstates+" in "+posteriorfiles[nfirstindex]);
		   }

	           while ((szLine = brposterior.readLine())!=null)
	           {
		      numlines++;
		   }
		   brposterior.close();	     
	    	
	           //loading in the posterior data
	           float[][] posterior = new float[numlines][numposteriorstates];
		
		   brposterior = Util.getBufferedReader(szposteriordir+"/"+szposteriorfiles_nfile);
		   int nline = 0;
	           brposterior.readLine();//specifies the header
	           brposterior.readLine();
	           while ((szLine = brposterior.readLine())!=null)
	           {
		      st = new StringTokenizer(szLine,"\t ");
		      float[] posterior_nline = posterior[nline];
                      for (int nstate = 0; nstate < numposteriorstates; nstate++)
                      {
		         float fval = Float.parseFloat(st.nextToken());
		         posterior_nline[nstate] = fval;
		         //keeps a tally of often each state occurs probabilistically
		         tallylabel[nstate] += fval;
		      }
		      nline++;
		   }
		   brposterior.close();

		   boolean[] counted = null;
 
		   BufferedReader brcoords = Util.getBufferedReader(szanchorpositions);
		   RecAnchorIndex theAnchorIndex = getAnchorIndex(szcolfields, busestrand, busesignal);
	           while ((szLine = brcoords.readLine())!=null)
	           {
		      if (szLine.trim().equals("")) continue;
	              //going through all the coordinates keeping only those for the current chromosome
	              String[] szLineA = szLine.split("\\s+");
                      String szreadchrom = szLineA[theAnchorIndex.nchromindex];        

	              if (szreadchrom.equals(szchrom))
	              {  
	                 int nanchor = (Integer.parseInt(szLineA[theAnchorIndex.npositionindex])-noffsetanchor);
	                 boolean bposstrand = true;
	                 if (busestrand)
	                 {
	                    String szstrand = szLineA[theAnchorIndex.nstrandindex];	    
	                    if (szstrand.equals("+"))
	                    {
			       bposstrand = true;
			    }
                            else if (szstrand.equals("-"))
                            {
			       bposstrand = false;
			    }
	                    else
	                    {
			       throw new IllegalArgumentException(szstrand +" is an valid strand. Strand should be '+' or '-'");
			    }     
			 }
		 
		         double damount;

	                 if ((busesignal)&&(theAnchorIndex.nsignalindex < szLineA.length))
	                 {
	                    damount = Double.parseDouble(szLineA[theAnchorIndex.nsignalindex]);
			 }
	                 else
                         {
	                    damount = 1;
			 }

	                 if (bposstrand)
	                 {
		            int ntallyindex = 0;
	                    for(int noffset= -numleft; noffset <= numright; noffset++)
	                    {
		               int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		               if ((nposindex >=0)&&(nposindex < posterior.length))
		               {
			          //within interval
		                  float[] posterior_nposindex = posterior[nposindex];
		                  double[] tallyoverlaplabel_ntallyindex = tallyoverlaplabel[ntallyindex];
		                  for (int nstate = 0; nstate < numposteriorstates; nstate++)
		                  {
				     //increments each state by the posterior amount
		                     tallyoverlaplabel_ntallyindex[nstate]+= damount*posterior_nposindex[nstate];		      
				  }
			       }
			       ntallyindex++;
			    }
			 }
	                 else
	                 {
		            int ntallyindex = 0;
	                    for(int noffset= numright; noffset >= -numleft; noffset--)
	                    {
		               int nposindex = (nanchor + nspacing*noffset)/nbinsize;

		               if ((nposindex >=0)&&(nposindex < posterior.length))
		               {
			          float[] posterior_nposindex = posterior[nposindex];
			          double[] tallyoverlaplabel_ntallyindex = tallyoverlaplabel[ntallyindex];
		                  for (int nstate = 0; nstate < numposteriorstates; nstate++)
		                  {
				     //increments each state by the posterior amount
	                             tallyoverlaplabel_ntallyindex[nstate]+=damount*posterior_nposindex[nstate];		      
				  }
			       }
			       ntallyindex++;		    			   
			    }
			 }
		      }
		   }		
		   brcoords.close();
		}
	    }
	}	

	if (!bposteriorfound)
	{
	    throw new IllegalArgumentException("No posterior file found for cell type "+szcell);
	}   
 
	outputneighborhood(tallyoverlaplabel,tallylabel,dsumoverlaplabel,szoutfile,nspacing,numright,
                           numleft,theColor,ChromHMM.convertCharOrderToStringOrder(chorder),sztitle,1,szlabelmapping,chorder);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * This handles outputing to text and heatmap from the fold enrichments tally info
     */
     private static void outputneighborhoodsignal(double[][] tallyoverlaplabel, double[] tallylabel, double[] dsumoverlaplabel, double dsumlabel, String szoutfile,
						  int nspacing, int numright, int numleft, Color theColor, String szstateorder,String sztitle,String szmarknames,
						  String szlabelmapping, char chorder) throws IOException
    {
        NumberFormat nf5 = NumberFormat.getInstance();
        nf5.setMaximumFractionDigits(5);
	nf5.setGroupingUsed(false);
	nf5.setMinimumFractionDigits(5);
 	    
	HashMap hmlabelExtend = makeLabelMapping(szlabelmapping);

	String[] collabels = new String[tallyoverlaplabel.length];
	System.out.println("Writing to file "+szoutfile+".txt");
	PrintWriter pw = new PrintWriter(new FileWriter(szoutfile+".txt"));

	//prints out the column labels
	//and stores them in collabels for the heatmap
	int ninterval = 0;
	pw.print("State ("+szstateorder+" order)");
	for(int npos = -nspacing*numleft; npos <= nspacing*numright; npos+=nspacing)
        {
	    pw.print("\t"+npos);
	    collabels[ninterval] = ""+npos;
	    ninterval++;
	}
	pw.println();

	//each row in the heatmap corresponds to a state and each column a position in tallyoverlaplabel
	double[][] heatmapfold = new double[tallyoverlaplabel[0].length][tallyoverlaplabel.length];
	int numelim = 0;

	//computes the fold enrichment for each non-empty state and outputs
	//the fold enrichment
	for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	{
	    if (tallylabel[nstate] > 0)
	    {
		//state actually occurs
		pw.print((nstate+1));
		String szsuffix;
		if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(nstate+1)))!=null)
		{
		    pw.print("_"+szsuffix);
		}

		for (int nfile = 0; nfile < tallyoverlaplabel.length; nfile++)
	        {
		    //the numerator is the fraction of the signal at this relative anchor assigned to this state
		    //the denominator is the fraction of the genome assigned to this state
	           double dfold = (tallyoverlaplabel[nfile][nstate]/dsumoverlaplabel[nfile])/
	                           (tallylabel[nstate]/dsumlabel);
		   heatmapfold[nstate][nfile] = dfold;
		   pw.print("\t"+nf5.format(dfold));
		}
		pw.println();
	    }
	    else
	    {
		//we are elminating this state when generating the heatmap since it was not found
		numelim++;
	    }
	}
	pw.close();

	StringTokenizer stheader = new StringTokenizer(szmarknames,"\t");
	String[] rowlabels = new String[tallyoverlaplabel[0].length];
	if (numelim > 0)
	{
	    //we need to collapse the heatmap
	   double[][] heatmapreduce = new double[heatmapfold.length-numelim][heatmapfold[0].length];
	   rowlabels = new String[heatmapreduce.length];
	   int nkeepindex = 0;
	   for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	   {
              if (tallylabel[nstate] > 0)
	      {
		  //keeping this index
		  for (int ncol = 0; ncol < heatmapfold[nstate].length; ncol++)
		  {
		      //copying over the contents
		      heatmapreduce[nkeepindex][ncol] = heatmapfold[nstate][ncol];
		  }
		  rowlabels[nkeepindex] = stheader.nextToken();//""+(nstate+1);
		  nkeepindex++;
	      }
	      stheader.nextToken();
	   }
	   heatmapfold = heatmapreduce;
	}
	else
	{
	    rowlabels = new String[tallyoverlaplabel[0].length];
	    for (int ni = 0; ni < rowlabels.length; ni++)
	    {
		rowlabels[ni] = ""+stheader.nextToken();
	    }
	}


        for (int ni = 0; ni < rowlabels.length; ni++)
	{
	   String szsuffix;
	   if ((szsuffix = (String) hmlabelExtend.get(""+chorder+rowlabels[ni]))!=null)
           {
	      rowlabels[ni] = rowlabels[ni]+"_"+szsuffix;
           }
        }

	makeEnrichmentHeatMap(heatmapfold, collabels, rowlabels,szoutfile,theColor,sztitle,"Position","State ("+szstateorder+" order)");
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * This handles outputing to text and heatmap form the fold enrichments tally info
     */
    private static void outputneighborhood(double[][] tallyoverlaplabel, double[] tallylabel, double[] dsumoverlaplabel, String szoutfile,
					   int nspacing, int numright, int numleft, Color theColor, String szstateorder,String sztitle,
                                           int noffset, String szlabelmapping, char chorder) throws IOException
    {
        NumberFormat nf5 = NumberFormat.getInstance();
        nf5.setMaximumFractionDigits(5);
	nf5.setGroupingUsed(false);
	nf5.setMinimumFractionDigits(5);
 	    
	HashMap hmlabelExtend = makeLabelMapping(szlabelmapping);

	//computes the total sum of signal each position, but summing over the 
	//signal at each state at that position
        for (int npos = 0; npos < tallyoverlaplabel.length; npos++)
        {
           for (int nindex = 0; nindex < tallyoverlaplabel[npos].length; nindex++)
           {
              dsumoverlaplabel[npos] += tallyoverlaplabel[npos][nindex];
           }
        }          

	double dsumlabel = 0;
	for (int ni = 0; ni < tallylabel.length; ni++)
        {
	    dsumlabel += tallylabel[ni];
	}

	String[] collabels = new String[tallyoverlaplabel.length];
	System.out.println("Writing to file "+szoutfile+".txt");
	PrintWriter pw = new PrintWriter(new FileWriter(szoutfile+".txt"));

	//prints out the column labels
	//and stores them in collabels for the heatmap
	int ninterval = 0;
	pw.print("State ("+szstateorder+" order)");
	for(int npos = -nspacing*numleft; npos <= nspacing*numright; npos+=nspacing)
        {
	    pw.print("\t"+npos);
	    collabels[ninterval] = ""+npos;
	    ninterval++;
	}
	pw.println();

	//each row in the heatmap corresponds to a state and each column a position in tallyoverlaplabel
	double[][] heatmapfold = new double[tallyoverlaplabel[0].length][tallyoverlaplabel.length];
	int numelim = 0;

	//computes the fold enrichment for each non-empty state and outputs
	//the fold enrichment
	for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	{
	    if (tallylabel[nstate] > 0)
	    {
		//state actually occurs
		pw.print((nstate+noffset));
                String szsuffix;
                if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(nstate+noffset)))!=null)
		{
	      	   pw.print("_"+szsuffix);
	        }

		for (int nfile = 0; nfile < tallyoverlaplabel.length; nfile++)
	        {
		    //the numerator is the fraction of the signal at this relative anchor assigned to this state
		    //the denominator is the fraction of the genome assigned to this state
	           double dfold = (tallyoverlaplabel[nfile][nstate]/dsumoverlaplabel[nfile])/
	                           (tallylabel[nstate]/dsumlabel);
		   heatmapfold[nstate][nfile] = dfold;
		   pw.print("\t"+nf5.format(dfold));
		}
		pw.println();
	    }
	    else
	    {
		//we are elminating this state when generating the heatmap since it was not found
		numelim++;
	    }
	}
	pw.close();


	String[] rowlabels = new String[tallyoverlaplabel[0].length];
	if (numelim > 0)
	{
	    //we need to collapse the heatmap
	   double[][] heatmapreduce = new double[heatmapfold.length-numelim][heatmapfold[0].length];
	   rowlabels = new String[heatmapreduce.length];
	   int nkeepindex = 0;
	   for (int nstate = 0; nstate < tallyoverlaplabel[0].length; nstate++)
	   {
              if (tallylabel[nstate] > 0)
	      {
		  //keeping this index
		  for (int ncol = 0; ncol < heatmapfold[nstate].length; ncol++)
		  {
		      //copying over the contents
		      heatmapreduce[nkeepindex][ncol] = heatmapfold[nstate][ncol];
		  }
		  rowlabels[nkeepindex] = ""+(nstate+noffset);
		  nkeepindex++;
	      }
	   }
	   heatmapfold = heatmapreduce;
	}
	else
	{
	    rowlabels = new String[tallyoverlaplabel[0].length];
	    for (int ni = 0; ni < rowlabels.length; ni++)
	    {
		rowlabels[ni] = ""+(ni+noffset);
	    }
	}

        for (int ni = 0; ni < rowlabels.length; ni++)
	{
     	   String szsuffix;
	   if ((szsuffix = (String) hmlabelExtend.get(chorder+rowlabels[ni]))!=null)
	   {
	      rowlabels[ni] = rowlabels[ni]+"_"+szsuffix;
           }
        }

	makeEnrichmentHeatMap(heatmapfold, collabels, rowlabels,szoutfile,theColor,sztitle,"Position","State ("+szstateorder+" order)");
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * szmainmodelfile is the reference set of emission parameters to use
     * szcomparedir is the directory for which the files beginning with emissions_ and ending with .txt will be compared
     * szoutputdir is the directory where the parametercorrelation .txt and .png output files will go
     * theColor is a Color object indicating the color of the heatmap
     */
    public static void makeModelEmissionCompare(String szmainmodelfile, String szcomparedir, String szoutputprefix,Color theColor) throws IOException
    {
	ArrayList alrecs = new ArrayList();

        NumberFormat nf5 = NumberFormat.getInstance();
        nf5.setMaximumFractionDigits(5);
	nf5.setGroupingUsed(false);
	nf5.setMinimumFractionDigits(5);

	File filemain = new File(szmainmodelfile);
	if ((!filemain.getName().startsWith("emissions_"))||(!szmainmodelfile.endsWith(".txt")))
        {
	    throw new IllegalArgumentException("Reference set of emission parameters must start with emissions_ and end with .txt");
	}
		
	BufferedReader bremissions = Util.getBufferedReader(szmainmodelfile);
	String szheader = bremissions.readLine();
	if (szheader == null)
	{
	    throw new IllegalArgumentException(szmainmodelfile+" is empty!");
	}
	StringTokenizer stheader = new StringTokenizer(szheader, "\t");
	int numcols = stheader.countTokens()-1;//substracts one since first column is state label column

	//maps each mark identifier to a consistent indext value
	HashMap hmNameToID = new HashMap();
	String szaxis = stheader.nextToken();
	int ncol = 0;
	while (stheader.hasMoreTokens())
	{
	    hmNameToID.put(stheader.nextToken(), Integer.valueOf(ncol));
	    ncol++;
	}

	//counts the number of states in the main file
	int numstatesmain = 0;
	String szLine;
	while ((szLine = bremissions.readLine())!=null)
        {
	   numstatesmain++;
	}
	bremissions.close();


	//this reads in the emission parameters of the reference model
	double[][] emissionparamsmain = new double[numstatesmain][numcols];
	bremissions = Util.getBufferedReader(szmainmodelfile);
	bremissions.readLine();
	String[] rowlabels = new String[numstatesmain];
	for (int nstate = 0; nstate <  numstatesmain; nstate++)
       	{
	    szLine = bremissions.readLine();
	    if (szLine == null)
	    {
		throw new IllegalArgumentException("Expecting "+numstatesmain+" lines in "+szmainmodelfile+" found fewer.");
	    }
	    StringTokenizer stLine = new StringTokenizer(szLine,"\t");
	    rowlabels[nstate] = stLine.nextToken();
	    for (ncol = 0; ncol < numcols; ncol++)
	    {
	       emissionparamsmain[nstate][ncol] = Double.parseDouble(stLine.nextToken());
	    }
	}
	bremissions.close();

	//this contains a mapping from mark identifier to regular column index
	int[] mappedcol = new int[hmNameToID.size()];

	File emissionsdir = new File(szcomparedir);
	String[] comparefiles = emissionsdir.list();

	for (int nfile = 0; nfile < comparefiles.length; nfile++)
	{
	    //compares to any file in the szcomparedir with the prefix "emissions_" and the suffixt ".txt"
	    RecEmissionFile theRecEmissionFile = new RecEmissionFile();
	    if ((comparefiles[nfile].startsWith("emissions_"))&&(comparefiles[nfile].endsWith(".txt")))
	    {
		//first figures out the number of states
		bremissions = Util.getBufferedReader(szcomparedir+"/"+comparefiles[nfile]);
		bremissions.readLine(); //remove header
		int numstates = 0;
		while ((szLine = bremissions.readLine())!=null)
		{
		    numstates++;
		}
		bremissions.close();
		theRecEmissionFile.numstates = numstates;
		theRecEmissionFile.szfilename = comparefiles[nfile];

		bremissions = Util.getBufferedReader(szcomparedir+"/"+comparefiles[nfile]);
		szheader = bremissions.readLine();
		if (szheader == null)
		{
		    throw new IllegalArgumentException(szcomparedir+"/"+comparefiles[nfile]+" is empty!");
		}
		stheader = new StringTokenizer(szheader, "\t");
	        int numcurrcols = stheader.countTokens()-1;
		stheader.nextToken();
		if (numcurrcols != mappedcol.length)
		{
		    throw new IllegalArgumentException("number of cols invalid found "+numcurrcols+" expecting "+mappedcol.length);
		}
		else
		{
		    //gets the mapping to the mark columns in this file to the reference order
	           ncol  = 0;
		   while (stheader.hasMoreTokens())
		   {
		       mappedcol[ncol] = ((Integer) hmNameToID.get(stheader.nextToken())).intValue();
		       ncol++;
		   }
		}

		theRecEmissionFile.emissionparams = new double[theRecEmissionFile.numstates][numcols];

		for (int nstate = 0; nstate <  theRecEmissionFile.numstates; nstate++)
		{
		    szLine = bremissions.readLine();
		    if (szLine == null)
		    {
			throw new IllegalArgumentException("Expecting "+theRecEmissionFile.numstates+" states in "+
							   szcomparedir+"/"+comparefiles[nfile]+" found fewer");
		    }
		    StringTokenizer stLine = new StringTokenizer(szLine,"\t");
		    stLine.nextToken(); //remove state ID
		    for (ncol = 0; ncol < numcols; ncol++)
		    {
			//storing the emission parameter for
			theRecEmissionFile.emissionparams[nstate][mappedcol[ncol]] = Double.parseDouble(stLine.nextToken());
		    }
		}
		//adds a record for this emission file containing the parameters, number of states, and number of marks
		alrecs.add(theRecEmissionFile);
		bremissions.close();
	    }
	}

	if (alrecs.size()==0)
	{
	    throw new IllegalArgumentException("No emission files to compare to found in "+szcomparedir);
	}
	RecEmissionFile[] theRecEmissionFileCompareA = new RecEmissionFile[alrecs.size()];
	double[][] bestcorrelation = new double[numstatesmain][alrecs.size()];

	//orders the comparison files by state and then ID
	for (int nindex = 0; nindex < theRecEmissionFileCompareA.length; nindex++)
	{
	    theRecEmissionFileCompareA[nindex] = (RecEmissionFile) alrecs.get(nindex);
	}
	Arrays.sort(theRecEmissionFileCompareA, new RecEmissionFileCompare());

	BufferedReader br;
	double dcorr;

	for (int nfile = 0; nfile < theRecEmissionFileCompareA.length; nfile++)
	{
	    //goes through each comparison file
	    int numstatescompare = theRecEmissionFileCompareA[nfile].numstates;
	    double[][] emissionparamscompare = theRecEmissionFileCompareA[nfile].emissionparams;
	    for (int nstatemain = 0; nstatemain < numstatesmain; nstatemain++)
	    {
		//goes through each main state
		double dbestcorr = -1;
		for (int nstatecompare = 0; nstatecompare < numstatescompare; nstatecompare++)
		{
		    //find the best correlated state in terms of emission parameters in the correlated vector
		   dcorr = Util.correlation(emissionparamsmain[nstatemain],emissionparamscompare[nstatecompare]);

		   if (dcorr > dbestcorr)
		   {
		       dbestcorr = dcorr;
		   }
		}
		//record this best correlation
		bestcorrelation[nstatemain][nfile] = dbestcorr;
	    }
	}


	String szfile = szoutputprefix+".txt";

	System.out.println("Writing to file "+szfile);
	PrintWriter pwcompare = new PrintWriter(new FileWriter(szfile));

	String[] collabels = new String[theRecEmissionFileCompareA.length];


	pwcompare.print("state");
	for (ncol = 0; ncol < theRecEmissionFileCompareA.length; ncol++)
	{
	    collabels[ncol] = ""+theRecEmissionFileCompareA[ncol].numstates;
	}

	for (ncol = 0; ncol < collabels.length; ncol++)
	{
	    pwcompare.print("\t"+theRecEmissionFileCompareA[ncol].szfilename);
	}
	pwcompare.println();

	//print to text file
	for (int nstate = 0; nstate < numstatesmain; nstate++)
	{
	   pwcompare.print(rowlabels[nstate]);
	   for (int nfile = 0; nfile < bestcorrelation[nstate].length; nfile++)
	   {
	       pwcompare.print("\t"+nf5.format(bestcorrelation[nstate][nfile]));
	   }
	   pwcompare.println();
	}
	pwcompare.close();


        HeatChart map = new HeatChart(bestcorrelation);

        map.setTitle("Best Emission Parameter Correlation");
        map.setXAxisLabel("Number of States in Model");
	map.setAxisValuesFont(new Font("SansSerif",0,20));
	map.setAxisLabelsFont(new Font("SansSerif",0,22));
	map.setTitleFont(new Font("SansSerif",0,24));
        map.setYAxisLabel(szaxis);
        map.setXValues(collabels);
        map.setYValues(rowlabels);
        map.setHighValueColour(theColor);
        map.saveToFile(new File(szoutputprefix+".png"));
        Util.printImageToSVG(map, szoutputprefix+".svg");
	System.out.println("Writing to file "+szoutputprefix+".png");
	System.out.println("Writing to file "+szoutputprefix+".svg");


    }
}

