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
 * This class handles generating the browser output.
 * The ChromHMM code was written by Jason Ernst.
 */

public class BrowserOutput
{


    /**
     * Record store an integer interval
     */
    static class BeginEndRec
    {
	int nbegin;
	int nend;
	BeginEndRec(int nbegin, int nend)
	{
	    this.nbegin = nbegin;
	    this.nend = nend;
	}
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *  Does a string comparison by first converting strings to integers
     *  If neither are integers does a stander string comparison
     */
    static class LabelCompare implements Comparator, Serializable
    {
	public int compare(Object o1, Object o2)
	{
	    String sz1 = (String) o1;
	    String sz2 = (String) o2;
            try
	    {
	       int n1 = Integer.parseInt(sz1);
	       int n2 = Integer.parseInt(sz2);
	       if (n1 < n2)
	       {
		   return -1;
	       }
	       else if (n1 > n2)
	       {
		   return 1;
	       }
	       else
	       {
		   return 0;
	       }
	    }
	    catch (NumberFormatException nfex)
	    {
		return (sz1.compareTo(sz2));
	    }
	}
    }

    /**
     * Stores the mapping from state IDs to R,G,B
     */
    HashMap hmcolor = new HashMap();

    /**
     * Stores the mapping from state IDs to descriptive labels
     */
    HashMap hmlabelExtend = new HashMap();

    /**
     * The file with the segmentation being loaded
     */
    String szsegmentfile;

    /**
     * Name of the file mapping color IDs to labels
     */
    String szcolormapping;

    /**
     * Name of the file mapping IDs to labels
     */
    String szidlabelmapping;

    /**
     * The name of the segmentation that will be used in the browser
     */
    String szsegmentationname;

    /**
     * The prefix for the output files
     */
    String szoutputfileprefix;

    /**
     * Number of states to make a segmentation for
     */
    int numstates;


    public BrowserOutput(String szsegmentfile, String szcolormapping,String szidlabelmapping, 
			 String szsegmentationname, String szoutputfileprefix, int numstates) throws IOException
    {
	this.szsegmentfile = szsegmentfile;
	this.szcolormapping =szcolormapping;
	this.szidlabelmapping = szidlabelmapping;
	this.szsegmentationname = szsegmentationname;
	this.szoutputfileprefix = szoutputfileprefix;
	this.numstates = numstates;

        hmcolor = new HashMap();
        hmlabelExtend = new HashMap();
	makeColorMapping();
	makeLabelMapping();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////


    private void makeColorMapping() throws IOException
    {
       String szLine;

       if (szcolormapping != null)
       {
          BufferedReader brcolor =  Util.getBufferedReader(szcolormapping);

          //Loading in a mapping from state ID to R,G,B color 

          while ((szLine = brcolor.readLine())!=null)
          {
             StringTokenizer st = new StringTokenizer(szLine,"\t");
	     String szID = st.nextToken();
	     String szColor = st.nextToken();
             hmcolor.put(szID, szColor);
	  }
	  brcolor.close();
       }
       else
       {
	   BufferedReader brsegment =  Util.getBufferedReader(szsegmentfile);
	   HashSet hs = new HashSet();
           while ((szLine =brsegment.readLine())!=null)
           {
              String[] szLineA = szLine.split("\\s+");
              hs.add(szLineA[3].substring(1)); //this removes first char giving ordering type
	   }
	   brsegment.close();

	   Iterator itr = hs.iterator();
	   int nmaxval=1;
	   int nminval=1;

           while (itr.hasNext())
	   {
	      int nval = Integer.parseInt((String)itr.next());

	      if (nval < nminval)
	      {
		  nminval = nval;
	      }
	      else if (nval > nmaxval)
	      {
		  nmaxval = nval;
	      }
	   }



	   nmaxval = Math.max(nmaxval, nminval + numstates - 1);
	   String[] szlabels = new String[nmaxval-nminval+1];

	   for (int nindex = nminval; nindex <= nmaxval; nindex++)
	   {
	      szlabels[nindex-nminval] = ""+nindex;
	   }
	   

	   int[][] colortable = getMaxDiverseOrder(szlabels.length);

	   for (int nindex = 0; nindex < szlabels.length; nindex++)
	   {
	       hmcolor.put(szlabels[nindex], colortable[nindex][0]+","+colortable[nindex][1]+","+colortable[nindex][2]);
	   }
       }
    }


    /**
     * From the set of 216 web colors excluding white and black selects the nsize most diverse
     * based on a greedy algorithm. From this set orders them by a greedy algorithm to find similiar
     * colors starting at Blue
     */
    private int[][] getMaxDiverseOrder(int nsize)
    {

	Random theRandom = new Random(132);
	int[] colors = new int[] {0, 51,102,153,204,255};
	int[][] colorscandidate =new int[colors.length*colors.length*colors.length][3];
	boolean[] taken = new boolean[colorscandidate.length];
	taken[0] = true;
	taken[taken.length-1] = true;

	//generating all pairs of r,g,b values from the above six
	int ncandidate = 0;
	for (int nr = 0; nr < colors.length; nr++)
	{
	    for (int ng = 0; ng < colors.length; ng++)
	    {
		for (int nb = 0; nb < colors.length; nb++)
		{
		    colorscandidate[ncandidate][0] = colors[nr];
		    colorscandidate[ncandidate][1] = colors[ng];
		    colorscandidate[ncandidate][2] = colors[nb];
		    ncandidate++;
		}
	    }
	}

	int[][] colorspaced = new int[nsize][3];
	colorspaced[0][0] = colorscandidate[colors.length-1][0];
	colorspaced[0][1] = colorscandidate[colors.length-1][1];
	colorspaced[0][2] = colorscandidate[colors.length-1][2];	

	for (int nstate =1; nstate < colorspaced.length; nstate++)
	{
	    int nbestdist = -1;
	    int nbestindex = -1;
	    for (ncandidate = 1; ncandidate < colorscandidate.length-1; ncandidate++)
	    {		
		//ignoring first and last color
		if (!taken[ncandidate])
		{
		    //not already used
		    int ndist;
		    int nmindist = Integer.MAX_VALUE;
		    //finding closest distance to already selected
		    for (int ncurrel = 0; ncurrel < nstate; ncurrel++)
		    {
		      
			ndist = Math.abs(colorscandidate[ncandidate][0]-colorspaced[ncurrel][0])+
                            Math.abs(colorscandidate[ncandidate][1]-colorspaced[ncurrel][1])+
			    Math.abs(colorscandidate[ncandidate][2]-colorspaced[ncurrel][2]);
 
		    
		       if (ndist < nmindist)
		       {
			  nmindist = ndist;
		       }
		    }
		
		    //if this candidate has the furthest closest distance so far then keep
		    if (nmindist > nbestdist)
		    {
		       nbestdist = nmindist;
		       nbestindex = ncandidate;
		    }
		}
	    }
	    if (nbestindex == -1)
	    {
		nbestindex = theRandom.nextInt(colorscandidate.length);
	    }
	    colorspaced[nstate][0] = colorscandidate[nbestindex][0];
	    colorspaced[nstate][1] = colorscandidate[nbestindex][1];
	    colorspaced[nstate][2] = colorscandidate[nbestindex][2];
	    taken[nbestindex] = true;
	}

	boolean[] takenordered = new boolean[colorspaced.length];
        int[][] colorsordered = new int[colorspaced.length][3];
	colorsordered[0][0] = colorspaced[0][0];
	colorsordered[0][1] = colorspaced[0][1];
	colorsordered[0][2] = colorspaced[0][2];
	//no ordering greedily taking the closest non selected one each iteration 
	for (int npos =1; npos < colorsordered.length; npos++)
	{
	    int nbestdist = Integer.MAX_VALUE;
	    int nbestindex = -1;
	    for (ncandidate = 1; ncandidate < colorsordered.length; ncandidate++)
	    {
		if (!takenordered[ncandidate])
		{
		   int ndist = Math.abs(colorspaced[ncandidate][0]-colorsordered[npos-1][0])+
		               Math.abs(colorspaced[ncandidate][1]-colorsordered[npos-1][1])+
		               Math.abs(colorspaced[ncandidate][2]-colorsordered[npos-1][2]);

		   if (ndist < nbestdist)
		   {
		      nbestdist = ndist;
		      nbestindex = ncandidate;
		   }
		}
	    }

	    takenordered[nbestindex] = true;
	    colorsordered[npos][0] = colorspaced[nbestindex][0];
	    colorsordered[npos][1] = colorspaced[nbestindex][1];
	    colorsordered[npos][2] = colorspaced[nbestindex][2];	    
	}

	return colorsordered;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////


    private void makeLabelMapping() throws IOException
    {
       if (szidlabelmapping != null)
       {
          BufferedReader bridlabel =  Util.getBufferedReader(szidlabelmapping);
          String szLine;

          //Loading in a mapping from state ID to  a label description

          while ((szLine = bridlabel.readLine())!=null)
          {
             StringTokenizer st = new StringTokenizer(szLine,"\t");
	     String szID = st.nextToken();
	     String szLabelExtend = st.nextToken();
	     hmlabelExtend.put(szID,szLabelExtend);
	  }	  
	  bridlabel.close();
       }

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Makes a single track browser view of the segmentation represented in szsegmentfile
     * szcolormapping is a two or three column text file which maps state ID to R,G,B color triples and optionally a state label
     * Name of segmentation in the browser file is given by szsegmentationame
     * Output is a file named szoutputfileprefix_browserdense.bed of segmentation viewable in a single track in dense mode in browser
     */
    public void makebrowserdense() throws IOException
    {
       System.out.println("Writing to file "+szoutputfileprefix+ChromHMM.SZBROWSERDENSEEXTENSION+".bed");
       PrintWriter pw = new PrintWriter(new FileWriter(szoutputfileprefix+ChromHMM.SZBROWSERDENSEEXTENSION+".bed"));

       BufferedReader brsegment =  Util.getBufferedReader(szsegmentfile);
       String szLine;
       boolean bfirst = true;

       while ((szLine =brsegment.readLine())!=null)
       {
	   StringTokenizer st = new StringTokenizer(szLine,"\t");
	   String szcurrchrom = st.nextToken();
	   int nbegin = Integer.parseInt(st.nextToken());
	   int nend = Integer.parseInt(st.nextToken());
	   String szFullID = st.nextToken();
	   String szID = szFullID.substring(1); //this removes ordering type
	   if (bfirst)
	   {
	       pw.println("track name=\""+szsegmentationname+"\" description=\" "+szsegmentationname+" ("+ChromHMM.convertCharOrderToStringOrder(szFullID.charAt(0))
                          +" ordered)"+"\" visibility=1 itemRgb=\"On\"");
	       bfirst = false;
	   }
	   String szColor = (String) hmcolor.get(szID);
 	   if (szColor == null)
           {
              throw new IllegalArgumentException("Color not given for "+szID);
           }

	   String szsuffix;
           if ((szsuffix = (String) hmlabelExtend.get(szFullID))!=null)
	   {
              szID = szID+"_"+szsuffix;
	   }

	   pw.println(szcurrchrom+"\t"+nbegin+"\t"+nend+"\t"+szID+"\t0\t.\t"+nbegin+"\t"+nend+"\t"+szColor);
        }
        brsegment.close();
        pw.close();      
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Makes a single track browser view of the segmentation represented in szsegmentfile
     * szcolormapping is a two or three column text file which maps state ID to R,G,B color triples and optionally a state label
     * Name of segmentation in the browser file is given by szsegmentationame
     * Output is a file named szoutputfileprefix_browserdense.bed of segmentation viewable with one state per row
     */
    public void makebrowserexpanded() throws IOException
    {
       System.out.println("Writing to file "+szoutputfileprefix+ChromHMM.SZBROWSEREXPANDEDEXTENSION+".bed");
       PrintWriter pw = new PrintWriter(new FileWriter(szoutputfileprefix+ChromHMM.SZBROWSEREXPANDEDEXTENSION+".bed"));

       String szLine;

       BufferedReader brsegment =  Util.getBufferedReader(szsegmentfile);

       //stores set of chromosomes and labels
       HashSet hschroms = new HashSet();
       HashSet hslabels = new HashSet();

       //stores for each chromosome the maximum coordinate
       HashMap hmchromMax = new HashMap();

       //stores the set of interval coordinates for each chromosome and label
       HashMap hmcoords = new HashMap();

       //maps a label without the prefix back to the full label
       HashMap hmlabelToFull = new HashMap();

       String szLabelFull=null;
       while ((szLine = brsegment.readLine())!=null)
       {
	   StringTokenizer st = new StringTokenizer(szLine,"\t");
	   String szchrom = st.nextToken();
	   int nbegin = Integer.parseInt(st.nextToken());
	   int nend = Integer.parseInt(st.nextToken());
	   szLabelFull = st.nextToken();
	   String szLabel = szLabelFull.substring(1);

	   hmlabelToFull.put(szLabel, szLabelFull);

	   hschroms.add(szchrom);
	   hslabels.add(szLabel);
	   ArrayList alRecs = (ArrayList) hmcoords.get(szchrom+"\t"+szLabel);
	   if (alRecs ==null)
	   {
	       //creating first entry for chromsome and coordinate
	       alRecs = new ArrayList();
	       hmcoords.put(szchrom+"\t"+szLabel,alRecs);
	   }
	   alRecs.add(new BeginEndRec(nbegin,nend));

	   //potentially updating maximum coordinate for chromosome
	   Object obj = ((Integer) hmchromMax.get(szchrom));	
	   if (obj != null)
           {
              int nval = ((Integer) obj).intValue();
	      
              if (nend > nval)
	      {
	         hmchromMax.put(szchrom,Integer.valueOf(nend));
              }
	   }
           else
           {  
              hmchromMax.put(szchrom,Integer.valueOf(nend));
	   }
       }
       brsegment.close();

      	
       //gets all the state labels and sorts them
       String[] szLabels = new String[hslabels.size()];
       Iterator itrLabels = hslabels.iterator();
       int nindex = 0;
       while (itrLabels.hasNext())
       {
	   szLabels[nindex] = (String) itrLabels.next();
	   nindex++;
       }
       Arrays.sort(szLabels, new LabelCompare());

       //gets all the chromsomes and sorts them
       String[] szChroms = new String[hschroms.size()];
       Iterator itrChroms = hschroms.iterator();
       nindex = 0;
       while (itrChroms.hasNext())
       {
	   szChroms[nindex] = (String) itrChroms.next();
	   nindex++;
       }
       Arrays.sort(szChroms);

       pw.println("track name=\"Expanded_"+szsegmentationname+"\" description=\" "+szsegmentationname+" ("+ChromHMM.convertCharOrderToStringOrder(szLabelFull.charAt(0))
                          +" ordered)"+"\" visibility=2 itemRgb=\"On\"");
       int nbrowserend = (int) (((Integer)hmchromMax.get(szChroms[0])).intValue()*.001)+1;
       pw.println("browser position "+szChroms[0]+":1-"+nbrowserend);

       for (int nlabel = szLabels.length-1; nlabel >=0; nlabel--)
       {
	   //UCSC browser seems to reverse the ordering of browser track files
	   String szcolor =  (String) hmcolor.get(""+szLabels[nlabel]);
	   for (int nchrom = 0; nchrom < szChroms.length; nchrom++)
	   {
	       //omits those segment labels not observed at all on chromosome
	       ArrayList alRecs  = (ArrayList) hmcoords.get(szChroms[nchrom]+"\t"+szLabels[nlabel]);
	       if (alRecs == null) continue;

               int nmax = ((Integer) hmchromMax.get(szChroms[nchrom])).intValue();

	       //this forces browser to display segment until the end of the chromosome
	       alRecs.add(new BeginEndRec(nmax-1,nmax));

	       int nsize = alRecs.size();
	       int nmin = ((BeginEndRec) alRecs.get(0)).nbegin;
	       int nfinalend = nmax;

	       String szoutlabel;
	       String szsuffix;
	       if ((szsuffix = (String) hmlabelExtend.get((String) hmlabelToFull.get(szLabels[nlabel])))!=null)
	       {
		   szoutlabel = szLabels[nlabel]+"_"+szsuffix;
	       }
	       else
	       {
		   szoutlabel = szLabels[nlabel];
	       }

	       pw.print(szChroms[nchrom]+"\t"+0+"\t"+nfinalend+"\t"+szoutlabel+"\t0\t.\t"+nmin+"\t"+nfinalend+"\t"+szcolor+"\t"+(nsize+1)+"\t");
	       pw.print(0); //forcing the display to start at the beginning of the chromosome
  	       for (int ni = 0; ni < nsize; ni++)
	       {
		   BeginEndRec theBeginEndRec = (BeginEndRec) alRecs.get(ni);
        	   int ndiff = theBeginEndRec.nend - theBeginEndRec.nbegin;
	           pw.print(",");
	           pw.print(ndiff);
	       }
	       pw.print("\t");
	       pw.print(0);
	       for (int ni = 0; ni < nsize; ni++)
               {
		   int nloc = ((BeginEndRec) alRecs.get(ni)).nbegin;
		   pw.print(",");
		   pw.print(nloc);
	       }
	       pw.println();
	   }
       }
       pw.close();
    }
}

