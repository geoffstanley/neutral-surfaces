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
package vgl.iisc.reebgraph.cmd;

import java.io.PrintStream;

import vgl.iisc.reebgraph.ui.BranchDecomposition;
import vgl.iisc.reebgraph.ui.ReebGraphData;
import vgl.iisc.reebgraph.ui.ReebGraphLayout;
import vgl.iisc.reebgraph.ui.ReebGraphLayout.Point;
import vgl.iisc.reebgraph.ui.SimplifyReebGraph;

public class LayoutReebGraph {

	ReebGraphData reebGraph;
	BranchDecomposition decomp;
	SimplifyReebGraph simplify;
	ReebGraphLayout layout;
	
	public void readReebGraph(String rgFile) {
		reebGraph = new ReebGraphData(rgFile);
		decomp = new BranchDecomposition(reebGraph, true);
		simplify = new SimplifyReebGraph(reebGraph, decomp);
		layout = new ReebGraphLayout();
	}
	
	public void writeLayoutCoordinates(String op, double sim) {
		simplify.simplify(sim);
		layout.layoutBranches(decomp, reebGraph, simplify.lastIncludedBranch);
		try {
			PrintStream pr = new PrintStream(op);
			int i = 0;
			for (Point node : layout.nodes) {
				if(simplify.drawNode[i]) {
					int v = reebGraph.nodes[i].v;
					pr.println(v + " " + node.loc.x + " " + node.loc.y + " " + node.loc.z);
				}
				i ++;
			}
			pr.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	
	public static void main(String[] args) {
		if(args.length != 2 && args.length != 3) {
			System.err.println("Invalid arguments");
			System.err.println("Arguments must be as follows : <input rg file> <output file> [simplification factor]");
			System.err.println("simplification factor is optional");
			System.exit(0);
		}
		String ip = args[0];
		String op = args[1];
		double sim = 0;
		if(args.length == 3) {
			try {
				sim = Double.parseDouble(args[2]);
				if(sim < 0 || sim > 1) {
					System.err.println("invalid simplification factor: should be a value between 0 and 1");
					System.exit(0);
				}
			} catch (Exception e) { 
				System.err.println("invalid simplification factor: should be a value between 0 and 1");
				System.exit(0);
			}
		}
		LayoutReebGraph layout = new LayoutReebGraph();
		layout.readReebGraph(ip);
		layout.writeLayoutCoordinates(op, sim);
	}	
}
