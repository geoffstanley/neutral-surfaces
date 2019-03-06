/*
 *	Copyright (C) 2012 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of Recon, a library to compute Reeb graphs.
 *
 *	Recon is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	Recon is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with Recon.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Author(s):	Harish Doraiswamy
 *	Version	 :	1.0
 *
 *	Modified by : -- 
 *	Date : --
 *	Changes  : --
 */
package vgl.iisc.cmd;


import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

import vgl.iisc.external.loader.DataLoader;
import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.recon.incore.*;
import vgl.iisc.reebgraph.ui.ReebGraphData;
import vgl.iisc.utils.Utilities;

public class ReCon {

	public static void main(String[] args) {
		
		try {
			Properties p = new Properties();
			p.load(new FileInputStream("input.properties"));
			String loaderType = p.getProperty("loader");
			String ip = p.getProperty("inputFile").trim();
			String fn = p.getProperty("inputFunction").trim();
			if(!Utilities.isInteger(fn)) {
				System.err.println("Input function should be a co-ordinate index (0 indicates given scalar function)");
				System.exit(1);
			}
			boolean adj = true;
			try {
				 adj = Boolean.parseBoolean(p.getProperty("adj").trim());	
			} catch (Exception e) {
				adj = true;
			}
			String op = null;
			try {
				 op = p.getProperty("output").trim();	
			} catch (Exception e) {
				op = null;
			}
			if(op != null && op.equalsIgnoreCase("")) {
				op = null;
			}
			String pFile = null;
			try {
				pFile = p.getProperty("partFile").trim();	
			} catch (Exception e) {
				pFile = null;
			}
			
			boolean aug = Boolean.parseBoolean(p.getProperty("augmentedRG").trim());
			MeshLoader loader = DataLoader.getLoader(loaderType);
			loader.setInputFile(ip);
			long st = System.nanoTime();
			
			if(aug) {
				ReconAlgorithmAug rg = new ReconAlgorithmAug();
				rg.useAdjacencies(adj);
				rg.computeReebGraph(loader, fn);
				if(op != null) {
					int []vmap = loader.getVertexMap();
					if(pFile != null) {
						rg.output(op, pFile);
						if(vmap!= null) {
							rewriteRG(op,vmap);
							rewritePart(pFile,vmap);
						}
					} else {
						rg.output(op);
						if(vmap!= null) {
							rewriteRG(op,vmap);
						}
					}
				}
			} else {
				ReConAlgorithmPrim rg = new ReConAlgorithmPrim();
				rg.useAdjacencies(adj);
				rg.computeReebGraph(loader, fn);
				if(op != null) {
					int []vmap = loader.getVertexMap();
					rg.output(op);
					if(vmap!= null) {
						rewriteRG(op,vmap);
					}
				}
			}
			long en = System.nanoTime();
			en -= st;
			double ms = en / 1000000;
			System.out.println("Total Time Taken : " + ms + " ms");
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	private static void rewritePart(String pFile, int[] vmap) {
		if(vmap == null) {
			return;
		}
		int [] invmap = new int[vmap.length];
		int valid = 0;
		for(int i = 0;i < vmap.length;i ++) {
			if(vmap[i] == -1) {
				continue;
			}
			invmap[vmap[i]] = i;
			valid = Math.max(valid, vmap[i]);
		}
		valid ++;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			int nv = Integer.parseInt(reader.readLine());
			int [] part = new int[nv];
			for(int i = 0;i < nv;i ++) {
				part[i] = Integer.parseInt(reader.readLine());
			}
			reader.close();
			PrintStream p = new PrintStream(pFile);
			p.println(valid);
			for(int i = 0;i < valid;i ++) {
				p.println(part[invmap[i]]);
			}
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void rewriteRG(String op, int[] vmap) {
		if(vmap == null) {
			return;
		}
		try {
			ReebGraphData rg = new ReebGraphData(op);
			PrintStream p = new PrintStream(op);
			p.println(rg.noNodes + " " + rg.noArcs);
			for(int i = 0;i < rg.noNodes;i ++) {
				int v = vmap[rg.nodes[i].v];
				p.println(v + " " + rg.nodes[i].fn + " " + rg.getTypeString(rg.nodes[i].type));
			}
			
			for(int i = 0;i < rg.noArcs;i ++) {
				p.print(rg.arcs[i].from + " " + rg.arcs[i].to + " ");
				p.println();
			}
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
