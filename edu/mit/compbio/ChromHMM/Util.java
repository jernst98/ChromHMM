
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
import java.util.zip.*;

import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;

import org.apache.batik.svggen.SVGGeneratorContext;
import org.tc33.jheatchart.HeatChart;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.dom.GenericDOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.DOMImplementation;



/**
 * Class contains general utility functions
 * The ChromHMM code was written by Jason Ernst 
 */
public class Util
{
    /**
     * Returns a buffered reader. If szFile ends in a ".gz" tries to open it as a gzip file
     * otherwise tries to open it as a normal file.
     */
    static BufferedReader getBufferedReader(String szFile) throws IOException
    {
	BufferedReader br;
       if (szFile.endsWith(".gz"))
       {
          br =new BufferedReader(new InputStreamReader(
			                  new GZIPInputStream(new FileInputStream(szFile))));
       }
       else
       {
          br = new BufferedReader(new FileReader(szFile));
       }
       return br;
    }

    /**
     * Computes the euclidean distance between the values in xvalues and the values in yvalues
     */
    static double euclid(double[] xvalues, double[] yvalues)
    {
       double ddist = 0;
       for (int ni =0; ni < xvalues.length; ni++)
       {
	   ddist += (xvalues[ni]-yvalues[ni])* (xvalues[ni]-yvalues[ni]);
       }

       return (Math.sqrt(ddist));
    }

    /**
     * Returns the correlation coeffiection between xvalues and yvalues
     * or 0 if either vector has variance 0
     */
    static double correlation(double[] xvalues, double[] yvalues)
    {
        double dsumx = 0,
	    dsumy = 0,
	    dsumxsq = 0,
	    dsumysq = 0,
	    dsumxy = 0,
	    dvarx,
	    dvary,
	    dcoeff;

        int numvalues = 0;

        for (int nindex = 0; nindex < xvalues.length; nindex++)
	{
	    dsumx += xvalues[nindex];
	    dsumy += yvalues[nindex];
	    dsumxsq += xvalues[nindex]*xvalues[nindex];
	    dsumysq += yvalues[nindex]*yvalues[nindex];
	    dsumxy  += xvalues[nindex]*yvalues[nindex];
	}
        numvalues = xvalues.length;

        if (numvalues==0)
        {
	    dcoeff = 0;
	}
        else
        {
	    dvarx = dsumxsq - dsumx*dsumx/numvalues;
	    dvary = dsumysq - dsumy*dsumy/numvalues;
	    double dvarxdvary = dvarx*dvary; 
	    if (dvarxdvary <= 0)
	    {
		dcoeff = 0;
	    }
	    else
	    {
       		dcoeff = (dsumxy - dsumx*dsumy/numvalues)/Math.sqrt(dvarxdvary);
	    }
	 }
         return dcoeff;
    }


    static void printImageToSVG(HeatChart map, String szoutfile) throws IOException
    {

        File svgFile = new File(szoutfile);
        Image imap = map.getChartImage();
        // Get a DOMImplementation.
        DOMImplementation domImpl =
	    GenericDOMImplementation.getDOMImplementation();
        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        // Create an instance of the SVG Generator.

	SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(document);
	ctx.setEmbeddedFontsOn(true);

        SVGGraphics2D svgGenerator = new SVGGraphics2D(ctx,true);
        //svgGenerator = canvas.getGraphics();

        svgGenerator.setColor(Color.white);

	map.getChartGraphics(svgGenerator);

        //svgGenerator.drawImage(imap, new AffineTransform(), null);
        //svgGenerator.drawImage(imap, new AffineTransform(1f,0f,0f,1f,0,0), null);
        //bimap1.paintComponent(svgGenerator);
        // Ask the test to render into the SVG Graphics2D implementation.
        //TestSVGGen test = new TestSVGGen();

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.
        boolean useCSS = true; // we want to use CSS style attributes
        Writer out = new OutputStreamWriter(new FileOutputStream(svgFile) , "UTF-8");

        svgGenerator.stream(out, useCSS);

    }
}
