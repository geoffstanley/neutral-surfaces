/*
 *	Copyright (C) 2010 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of libRG, a library to compute Reeb graphs.
 *
 *	Recon is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	libRG is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with libRG.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Author(s):	Harish Doraiswamy
 *	Version	 :	1.0
 *
 *	Modified by : Geoffrey J. Stanley 
 *	Date : 25/11/2018
 *	Changes  : Commented out "No. of Nodes" println
 */
package vgl.iisc.recon.incore;

import java.io.PrintStream;
import java.io.Serializable;


public class StoreReebGraph implements Serializable {

	private static final long serialVersionUID = 1L;

	int [] nodes;
	byte [] nodeType;
	double [] nodeFn;
	MyIntList xarc;
	MyIntList yarc;

	public StoreReebGraph(int no) {
		nodes = new int[no];
		nodeType = new byte[no];
		nodeFn = new double[no];
		xarc = new MyIntList((int) (no*1.5));
		yarc = new MyIntList((int) (no*1.5));
	}
	
	int curNode = 0;
	public void addNode(int v, double fn, byte type) {
		nodes[curNode] = v;
		nodeType[curNode] = type;
		nodeFn[curNode] = fn;
		curNode ++;
	}
	
	CleanReebGraph rg;
	public void setup() {
		 rg = new CleanReebGraph();
		 // System.out.println("No. of Nodes: " + curNode); // GJS commented out
		 for(int i = 0;i < curNode;i ++) {
			 rg.addNode(nodes[i], nodeFn[i], nodeType[i]);
		 }
		 
		 for(int i = 0;i < xarc.length;i ++) {
			 rg.addArc(xarc.array[i], yarc.array[i]);
		 }
		 rg.setup();
	}
	
	public void addArc(int v1, int v2) {
		xarc.add(v1);
		yarc.add(v2);
	}

	public void removeDeg2Nodes() {
		rg.removeDeg2Nodes();
	}
	
	
	public void outputReebGraph(PrintStream p) {
		rg.outputReebGraph(p);
	}
	
	public void outputReebGraph(PrintStream p, PrintStream part, int noVerts) {
		rg.outputReebGraph(p);
		rg.outputSegmentation(part, noVerts);
	}
}
