/**
 * ChromHMM - automating chromatin state discovery and characterization
 * Copyright (C) 2008-2012 Massachusetts Institute of Technology
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

package edu.mit.compbio.ChromHMM;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

public class ConvertGeneTable
{

    /**
     * Reads in a RefSeq or equivalently formatted gene table from the UCSC browser and outputs converted files
     * the columns in order from the RefSeq gene table are (additional subsequent columns are ignored if present)
     * bin
     * name
     * chrom
     * strand
     * txStart
     * txEnd
     * cdsStart
     * cdsEnd
     * exonCount
     * exonStarts
     * exonEnds
     

     * or a bigGenePred format, where only the first 12 columns are considered:
     * chrom
     * chromStart
     * chromEnd
     * name
     * score
     * strand
     * thickStart
     * thickEnd
     * reserved
     * blockCount
     * blockSizes
     * chromStarts
    */

    static void convertGeneTableToAnnotations(String sztable, String szprefix,
                                              String szassembly, String szcoorddir, 
					      String szanchordir, String szchromlengths, 
                                              int npromoterwindow, boolean bgzip, 
                                              boolean bnobin, boolean bnoheader, boolean bbiggenepred) throws IOException
    {

	//String sztable = args[0];
	//String szassembly = args[1];
	//String szoutdir = args[2];
	//String szanchordir = args[3];
	//String szchromlengths = args[4];
	String szLine;
	HashMap hmlengths = new HashMap();
	//BufferedReader brlength = new BufferedReader(new FileReader(szchromlengths));
	BufferedReader brlength = Util.getBufferedReader(szchromlengths);
	while ((szLine = brlength.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t ");
            if (st.countTokens()==0)
	    {
	       throw new IllegalArgumentException("empty line found in "+szchromlengths);
	    }
            else if (st.countTokens()==1)
	    {
	       throw new IllegalArgumentException("Only one token found in line "+szLine+" of "+szchromlengths);
	    }
	    hmlengths.put(st.nextToken(),Integer.valueOf(st.nextToken()));
	}
	brlength.close();
        BufferedReader br = Util.getBufferedReader(sztable);


	PrintWriter pwtss = null;
	PrintWriter pwanchortss = null;
	PrintWriter pwanchortes = null;
	PrintWriter pwgene = null;
	PrintWriter pwexon = null;
	PrintWriter pwtes = null;
	PrintWriter pwtss2kb = null;

	GZIPOutputStream pwtsszip = null;
        GZIPOutputStream pwanchortsszip = null;
        GZIPOutputStream pwanchorteszip = null;
        GZIPOutputStream pwgenezip = null;
        GZIPOutputStream pwexonzip = null;
        GZIPOutputStream pwteszip = null;
        GZIPOutputStream pwtss2kbzip = null;

	if (bgzip)
	{
	   pwtsszip =  new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"TSS."+szassembly+".bed.gz"));
	   pwanchortsszip = new GZIPOutputStream(new FileOutputStream(szanchordir+"/"+szprefix+"TSS."+szassembly+".txt.gz"));
	   pwanchorteszip = new GZIPOutputStream(new FileOutputStream(szanchordir+"/"+szprefix+"TES."+szassembly+".txt.gz"));
	   pwgenezip = new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"Gene."+szassembly+".bed.gz"));

	   pwexonzip = new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"Exon."+szassembly+".bed.gz"));
	   pwteszip = new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"TES."+szassembly+".bed.gz"));

	   if (npromoterwindow % 1000 == 0)
	   {
	       pwtss2kbzip = new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"TSS"+(npromoterwindow/1000)+"kb."+szassembly+".bed.gz"));
	   }
	   else
	   {
	       pwtss2kbzip = new GZIPOutputStream(new FileOutputStream(szcoorddir+"/"+szprefix+"TSS"+npromoterwindow+"bp."+szassembly+".bed.gz"));
	   }
	}
	else
	{
	   pwtss = new PrintWriter(szcoorddir+"/"+szprefix+"TSS."+szassembly+".bed");
	   pwanchortss = new PrintWriter(szanchordir+"/"+szprefix+"TSS."+szassembly+".txt");
	   pwanchortes = new PrintWriter(szanchordir+"/"+szprefix+"TES."+szassembly+".txt");
	   pwgene = new PrintWriter(szcoorddir+"/"+szprefix+"Gene."+szassembly+".bed");

	   pwexon = new PrintWriter(szcoorddir+"/"+szprefix+"Exon."+szassembly+".bed");
	   pwtes = new PrintWriter(szcoorddir+"/"+szprefix+"TES."+szassembly+".bed");
	 
	   if (npromoterwindow % 1000 == 0)
	   {
	      pwtss2kb = new PrintWriter(szcoorddir+"/"+szprefix+"TSS"+(npromoterwindow/1000)+"kb."+szassembly+".bed");
	   }
	   else
	   {
	      pwtss2kb = new PrintWriter(szcoorddir+"/"+szprefix+"TSS"+npromoterwindow+"bp."+szassembly+".bed");
	   }
	}

	HashSet hstss = new HashSet();
	HashSet hsanchortss = new HashSet();
	HashSet hsanchortes = new HashSet();
	HashSet hsgene = new HashSet();
	HashSet hsexon = new HashSet();
	HashSet hstes = new HashSet();
	HashSet hstss2kb = new HashSet();

	if (!bnoheader)
	{
	   br.readLine();
	}

	while ((szLine = br.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t",true);
	    String szchrom = null;
	    String szstrand = null;
	    String sztxStart = null;
	    String sztxEnd = null;
	    String szexonStarts = null;
	    String szexonEnds = null;
	    String szexonSizes = null;
	    String szchromStarts = null;

	    if (bbiggenepred)
	    {
	       szchrom = st.nextToken().trim();
	       if (!szchrom.equals("\t"))
	          st.nextToken().trim();
	       sztxStart = st.nextToken().trim();
	       if (!sztxStart.equals("\t"))
	          st.nextToken().trim();
	       sztxEnd = st.nextToken().trim();
	       if (!sztxEnd.equals("\t"))
	          st.nextToken().trim();
	       String szname = st.nextToken().trim();
	       if (!szname.equals("\t"))
	          st.nextToken().trim();
	       String szscore = st.nextToken().trim();
	       if (!szscore.equals("\t"))
	          st.nextToken().trim();
	       szstrand = st.nextToken().trim();
	       if (!szstrand.equals("\t"))
	          st.nextToken().trim();
	       String szthickStart = st.nextToken().trim();
	       if (!szthickStart.equals("\t"))
	          st.nextToken().trim();
	       String szthickEnd = st.nextToken().trim();
	       if (!szthickEnd.equals("\t"))
	          st.nextToken().trim();
	       String szreserved = st.nextToken().trim();
	       if (!szreserved.equals("\t"))
	          st.nextToken().trim();
	       String szexonCount = st.nextToken().trim();
	       if (!szexonCount.equals("\t"))
	          st.nextToken().trim();
	       szexonSizes = st.nextToken().trim();
	       if (!szexonSizes.equals("\t"))
	          st.nextToken().trim();
	       szchromStarts = st.nextToken().trim();
	    }
	    else
	    {
	       if (!bnobin)
	       {
		   //updated in 1.25
		  String szbin = null;
	          szbin = st.nextToken().trim();
	          if (!szbin.equals("\t"))
	             st.nextToken().trim();
	       }
	       String szname = st.nextToken().trim();
	       if (!szname.equals("\t"))
	          st.nextToken().trim();
	       szchrom = st.nextToken().trim();
	       if (!szchrom.equals("\t"))
	          st.nextToken().trim();
	       szstrand = st.nextToken().trim();
	       if (!szstrand.equals("\t"))
	          st.nextToken().trim();
	       sztxStart = st.nextToken().trim();
	       if (!sztxStart.equals("\t"))
	          st.nextToken().trim();
	       sztxEnd = st.nextToken().trim();
	       if (!sztxEnd.equals("\t"))
	          st.nextToken().trim();
	       String szcdsStart = st.nextToken().trim();
	       if (!szcdsStart.equals("\t"))
	          st.nextToken().trim();
	       String szcdsEnd = st.nextToken().trim();
	       if (!szcdsEnd.equals("\t"))
	          st.nextToken().trim();
	       String szexonCount = st.nextToken().trim();
	       if (!szexonCount.equals("\t"))
	          st.nextToken().trim();
	       szexonStarts = st.nextToken().trim();
	       if (!szexonStarts.equals("\t"))
	          st.nextToken().trim();
	       szexonEnds = st.nextToken().trim();

	    }
	    //removed 1.22
	    //if (!szexonEnds.equals("\t"))
	    //	st.nextToken();
	    
	    //String szscore = st.nextToken();
	    //if (!szscore.equals("\t"))
	    //st.nextToken();
	    //String szname2 = st.nextToken();
	    //if (!szname2.equals("\t"))
	    //st.nextToken();
	    //String szcdsStartStat = st.nextToken();
	    //if (!szcdsStartStat.equals("\t"))
	    //	st.nextToken();
	    //String szcdsEndStat = st.nextToken();
	    //if (!szcdsEndStat.equals("\t"))
	    //st.nextToken();
	    //String szexonFrames = st.nextToken();


	    int ntes=-1;
	    int ntss=-1;
	    if (szstrand.equals("+"))
	    {
	       ntss = Integer.parseInt(sztxStart);
	       ntes = Integer.parseInt(sztxEnd)-1;
	    }
	    else if (szstrand.equals("-"))
	    {
	       ntss = Integer.parseInt(sztxEnd)-1;
	       ntes = Integer.parseInt(sztxStart);
	    }
	    else
	    {
		throw new IllegalArgumentException("invalid strand\t"+szstrand);
	    }

	    //updated in 1.18 to give better error message if chromosome not found
	    Integer objchromlength = ((Integer) hmlengths.get(szchrom));
	    int nchromlength = -1;

	    if (objchromlength != null)
	    {
		nchromlength = objchromlength.intValue();
	    }
	    else
	    {
                throw new IllegalArgumentException("did not find length in chromosome file for chromosome: "+szchrom);
		// =((Integer) hmlengths.get(szchrom)).intValue();
	    }

	    String sztssOut = szchrom+"\t"+ntss+"\t"+(ntss+1);
	    String sztesOut = szchrom+"\t"+ntes+"\t"+(ntes+1);
	    String sztss2kbOut = szchrom+"\t"+Math.max((ntss-npromoterwindow),0)+"\t"+Math.min(ntss+npromoterwindow+1,nchromlength);
	    String szgeneOut = szchrom+"\t"+sztxStart+"\t"+sztxEnd;
	    String szanchorTSSOut = szchrom+"\t"+ntss+"\t"+szstrand;
	    String szanchorTESOut = szchrom+"\t"+ntes+"\t"+szstrand;
	    //added two allow " quotes for exon list
	    StringTokenizer stexonStarts = null;
	    StringTokenizer stexonEnds = null;
	    StringTokenizer stexonSizes = null;
	    StringTokenizer stchromStarts = null;

	    if (bbiggenepred)
	    {
	       stexonSizes = new StringTokenizer(szexonSizes,",\"");
	       stchromStarts = new StringTokenizer(szchromStarts,",\"");
	    }
	    else
	    {
	       stexonStarts = new StringTokenizer(szexonStarts,",\"");
	       stexonEnds = new StringTokenizer(szexonEnds,",\"");
	    }

	    if (bgzip)
	    {
               if (!hstss.contains(sztssOut))
	       {
		  byte[] btformat = (sztssOut+"\n").getBytes();
		  pwtsszip.write(btformat,0,btformat.length);
		  //pwtss.println(sztssOut);
		  hstss.add(sztssOut);	       
	       }

	       if (!hsanchortss.contains(szanchorTSSOut))
	       {
		   byte[] btformat = (szanchorTSSOut+"\n").getBytes();
		  pwanchortsszip.write(btformat,0,btformat.length);

		  //pwanchortss.println(szanchorTSSOut);
		  hsanchortss.add(szanchorTSSOut);
	       }

	       if (!hsanchortes.contains(szanchorTESOut))
	       {
		   byte[] btformat = (szanchorTESOut+"\n").getBytes();
		  pwanchorteszip.write(btformat,0,btformat.length);

		  //pwanchortes.println(szanchorTESOut);
		  hsanchortes.add(szanchorTESOut);
	       }

	       if (!hstes.contains(sztesOut))
	       {
		   byte[] btformat = (sztesOut+"\n").getBytes();
		  pwteszip.write(btformat,0,btformat.length);

	          //pwtes.println(sztesOut);
	          hstes.add(sztesOut);
	       }

	       if (!hstss2kb.contains(sztss2kbOut))
	       {
		   byte[] btformat = (sztss2kbOut+"\n").getBytes();
		  pwtss2kbzip.write(btformat,0,btformat.length);

	          //pwtss2kb.println(sztss2kbOut);
	          hstss2kb.add(sztss2kbOut);
	       }

	       if (!hsgene.contains(szgeneOut))
	       {
		   byte[] btformat = (szgeneOut+"\n").getBytes();
		  pwgenezip.write(btformat,0,btformat.length);

	          //pwgene.println(szgeneOut);
	          hsgene.add(szgeneOut);
	       }

	       StringTokenizer stStarts = null;
	       if (bbiggenepred)
	       {
	          stStarts = stchromStarts;
	       }
	       else
	       {
	          stStarts = stexonStarts;
	       }

	       while (stStarts.hasMoreTokens())
	       {
		  String szexonOut = null;
		  if (bbiggenepred)
		  {
		     int nchromstart = Integer.parseInt(sztxStart);
		     int nexonchromstart = Integer.parseInt(stStarts.nextToken().trim());
		     int nexonsize = Integer.parseInt(stexonSizes.nextToken().trim());
		     szexonOut = szchrom+"\t"+(nchromstart+nexonchromstart)+"\t"+(nchromstart+nexonchromstart+nexonsize)+"\n";
		  }
		  else
		  {
		     szexonOut = szchrom+"\t"+stStarts.nextToken().trim()+"\t"+stexonEnds.nextToken().trim()+"\n";
		  }

		  if (!hsexon.contains(szexonOut))
		  {
                     byte[] btformat = szexonOut.getBytes();
                     pwexonzip.write(btformat,0,btformat.length);

		     hsexon.add(szexonOut);
		  }

		  //pwexon.println(szexonOut);
	       }
	    }
	    else
	    {
               if (!hstss.contains(sztssOut))
	       {
		  pwtss.println(sztssOut);
		  hstss.add(sztssOut);	       
	       }

	       if (!hsanchortss.contains(szanchorTSSOut))
	       {
		  pwanchortss.println(szanchorTSSOut);
		  hsanchortss.add(szanchorTSSOut);
	       }

	       if (!hsanchortes.contains(szanchorTESOut))
	       {
		  pwanchortes.println(szanchorTESOut);
		  hsanchortes.add(szanchorTESOut);
	       }

	       if (!hstes.contains(sztesOut))
	       {
	          pwtes.println(sztesOut);
	          hstes.add(sztesOut);
	       }

	       if (!hstss2kb.contains(sztss2kbOut))
	       {
	          pwtss2kb.println(sztss2kbOut);
	          hstss2kb.add(sztss2kbOut);
	       }

	       if (!hsgene.contains(szgeneOut))
	       {
	          pwgene.println(szgeneOut);
	          hsgene.add(szgeneOut);
	       }

	       StringTokenizer stStarts = null;
	       if (bbiggenepred)
	       {
	          stStarts = stchromStarts;
	       }
	       else
	       {
	          stStarts = stexonStarts;
	       }

	       while (stStarts.hasMoreTokens())
	       {

		  String szexonOut = null;
		  if (bbiggenepred)
		  {
		     int nchromstart = Integer.parseInt(sztxStart);
		     int nexonchromstart = Integer.parseInt(stStarts.nextToken().trim());
		     int nexonsize = Integer.parseInt(stexonSizes.nextToken().trim());
		     szexonOut = szchrom+"\t"+(nchromstart+nexonchromstart)+"\t"+(nchromstart+nexonchromstart+nexonsize);
		  }
		  else
		  {
		     szexonOut = szchrom+"\t"+stStarts.nextToken().trim()+"\t"+stexonEnds.nextToken().trim();
		  }

                  if (!hsexon.contains(szexonOut))
		  {
		     pwexon.println(szexonOut);
		     hsexon.add(szexonOut);
		  }
	       }
	    }
	}
	br.close();

	if (bgzip)
	{
	   pwtsszip.finish(); //tss
	   pwgenezip.finish(); //gene
	   pwanchortsszip.finish();
	   pwanchorteszip.finish();
	   pwexonzip.finish(); //exon
	   pwteszip.finish(); //tss
	   pwtss2kbzip.finish(); //2kb


	   pwtsszip.close(); //tss
	   pwgenezip.close(); //gene
	   pwanchortsszip.close();
	   pwanchorteszip.close();
	   pwexonzip.close(); //exon
	   pwteszip.close(); //tss
	   pwtss2kbzip.close(); //2kb
	}
	else
	{
	   pwtss.close(); //tss
	   pwgene.close(); //gene
	   pwanchortss.close();
	   pwanchortes.close();
	   pwexon.close(); //exon
	   pwtes.close(); //tss
	   pwtss2kb.close(); //2kb
	}

    }

}
