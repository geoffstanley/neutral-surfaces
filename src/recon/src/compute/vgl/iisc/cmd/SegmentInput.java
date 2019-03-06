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
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;

import vgl.iisc.cmd.WritePartitionedOutput;
import vgl.iisc.reebgraph.ui.ReebGraphData;

public class SegmentInput {

	ReebGraphData rgData;
	public int [] nodeComp;
	
	public void getPartitionInput(String partFile, String rgFile) throws IOException {
		rgData = new ReebGraphData(rgFile);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(partFile));
			int nv = Integer.parseInt(reader.readLine().trim());
			nodeComp = new int[nv];
			for(int i = 0;i < nv;i ++) {
				nodeComp[i] = Integer.parseInt(reader.readLine().trim());
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void output(String op) {
		try {
			ObjectOutputStream p = new ObjectOutputStream(new FileOutputStream(op));
			p.writeObject(nodeComp);
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void outputVRF(String op, String tetFile, String fn) {
		WritePartitionedOutput p = new WritePartitionedOutput();
		p.writePartition(tetFile, fn, rgData, nodeComp, op);
	}
	
	public void outputVRFFromRaw(String op, String rawFile, boolean ascii) {
		WritePartitionedOutput p = new WritePartitionedOutput();
		p.setRaw(true);
		p.writePartition(rawFile, nodeComp, op, ascii);
	}
	
	public void outputOFF(String offFile, String ip) {
		System.gc();
		WritePartitionedOutput p = new WritePartitionedOutput();
		p.writeOffFile(ip, offFile, nodeComp);
	}
}
