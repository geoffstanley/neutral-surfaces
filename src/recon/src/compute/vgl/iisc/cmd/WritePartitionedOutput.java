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

import static vgl.iisc.utils.Utilities.er;
import static vgl.iisc.utils.Utilities.pr;
import static vgl.iisc.utils.Utilities.splitString;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import vgl.iisc.cmd.data.GridReader;
import vgl.iisc.cmd.data.TriangleDataOpt;
import vgl.iisc.external.loader.DataLoader;
import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;
import vgl.iisc.reebgraph.ui.ReebGraphData;
import vgl.iisc.reebgraph.ui.ReebGraphData.Arc;
import vgl.iisc.utils.Utilities;

public class WritePartitionedOutput {
	class Edges {
		ArrayList<Integer> edges = new ArrayList<Integer>();
	}
	
	class Tetrahedron {
		int [] v = new int[4];
	}
	
	Tetrahedron [] tets;
	TriangleDataOpt data;
	ReebGraphData rgData;
	int [] arcs;
	static short size = 64;

	
	public void writePartition(String tetFile, String fn, String rgFile, String opFile, String partFile) {
		MeshLoader loader = DataLoader.getLoader("TET");
		loader.setInputFile(tetFile);
		loadData(loader, fn);
		rgData = new ReebGraphData(rgFile);
		readPartition(partFile);
		
		writePartition(opFile);
	}
	
	public void writePartition(String tetFile, String fn, ReebGraphData rgData, int [] nodeComp, String opFile) {
		MeshLoader loader = DataLoader.getLoader("TET");
		loader.setInputFile(tetFile);
		loadData(loader, fn);
		this.rgData = rgData;
		arcs = nodeComp;
		
		writePartition(opFile);
	}

	private void writePartition(String opFile) {
		sampleTets();
		writeVRFFile(opFile);
	}
	
	private void readPartition(String partFile) {
		try {
			ObjectInputStream p = new ObjectInputStream(new FileInputStream(partFile));
			arcs = (int[]) p.readObject();
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void loadData(MeshLoader loader, String ftype) {
		data = new TriangleDataOpt(true);
		ArrayList<Tetrahedron> teta = new ArrayList<Tetrahedron>(); 
		try {
			int nv = loader.getVertexCount();
			data.setNoOfVertices(nv);
			Simplex sim = loader.getNextSimplex();
			while (sim != null) {
				if (sim instanceof vgl.iisc.external.types.Vertex) {
					vgl.iisc.external.types.Vertex v = (vgl.iisc.external.types.Vertex) sim;
					double fn;
					int f = Integer.parseInt(ftype);
					if (f == 0) {
						fn = v.f;
					} else {
						fn = v.c[f - 1];
					}
					if(v.c.length < 3) {
						data.addVertex(v.c[0], v.c[1], 0, fn);	
					} else {
						data.addVertex(v.c[0], v.c[1], v.c[2], fn);
					}
					
				} else if (sim instanceof vgl.iisc.external.types.Edge) {
					pr("Simplicial meshes with edges are not yet supported. This edge will be ignored");
				} else if (sim instanceof vgl.iisc.external.types.Triangle) {
					vgl.iisc.external.types.Triangle t = (vgl.iisc.external.types.Triangle) sim;
					data.addTriangle(t.v1, t.v2, t.v3);
				} else if (sim instanceof vgl.iisc.external.types.Tetrahedron) {
					vgl.iisc.external.types.Tetrahedron t = (vgl.iisc.external.types.Tetrahedron) sim;
					data.addTetraHedra(t.v1, t.v2, t.v3, t.v4);
					Tetrahedron tt = new Tetrahedron();
					tt.v[0] = t.v1;
					tt.v[1] = t.v2;
					tt.v[2] = t.v3;
					tt.v[3] = t.v4;
					teta.add(tt);
				} else {
					er("Invalid Simplex");
				}
				sim = loader.getNextSimplex();
			}
			pr("Finished reading data from file. Loading it......");
			data.setupTriangles();
			pr("Successfully loaded Data");
			tets = teta.toArray(new Tetrahedron[0]);
			normalizeVertices();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	private void normalizeVertices() {
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		
		double minx = Double.MAX_VALUE;
		double miny = Double.MAX_VALUE;
		double minz = Double.MAX_VALUE;
		
		double maxx = - Double.MAX_VALUE;
		double maxy = - Double.MAX_VALUE;
		double maxz = - Double.MAX_VALUE;
		
		for(int i = 0;i < data.coords.length;i ++) {
			min = Math.min(min, data.fnVertices[i]);
			max = Math.max(max, data.fnVertices[i]);
			
			minx = Math.min(minx, data.coords[i][0]);
			miny = Math.min(miny, data.coords[i][1]);
			minz = Math.min(minz, data.coords[i][2]);
			
			maxx = Math.max(maxx, data.coords[i][0]);
			maxy = Math.max(maxy, data.coords[i][1]);
			maxz = Math.max(maxz, data.coords[i][2]);
		}
		
		double extentx = maxx - minx;
		double extenty = maxy - miny;
		double extentz = maxz - minz;
		
		
		double extent = Math.max(extentx,extenty);
		extent = Math.max(extent,extentz);
		double scale = size / extent;
		
		double diff = max - min;
		for(int i = 0;i < data.coords.length;i ++) {
			if(diff != 0) {
				data.fnVertices[i] = (data.fnVertices[i] - min) / diff;
			} else {
				data.fnVertices[i] = 1;
			}
			
			data.coords[i][0] -= minx;
			data.coords[i][1] -= miny;
			data.coords[i][2] -= minz;
			
			data.coords[i][0] += (extent - extentx)/2;
			data.coords[i][1] += (extent - extenty)/2;
			data.coords[i][2] += (extent - extentz)/2;
			
			data.coords[i][0] *= scale;
			data.coords[i][1] *= scale;
			data.coords[i][2] *= scale;
		}
	}

	double [][][] gData;
	
	private void sampleTets() {
		System.out.println("Converting tets to grid");
		
		int [][][] ct = new int[size][size][size];
		int [][][] cl = new int[size][size][size];
		gData = new double [size][size][size];
		double max = -1;
		int nt = tets.length;
		int xx = nt / 100;
		for(int i = 0;i < nt;i ++) {
			if(i % xx == 0) {
				System.out.println(i / xx);
			}
			
			// sample Tet
			
			double minx;
			double miny;
			double minz;
			
			double maxx;
			double maxy;
			double maxz;
			
			minx = maxx = data.coords[tets[i].v[0]][0];
			miny = maxy = data.coords[tets[i].v[0]][1];
			minz = maxz = data.coords[tets[i].v[0]][2];
			
			for(int j = 1;j < 4;j ++) {
				minx = Math.min(minx, data.coords[tets[i].v[j]][0]);
				miny = Math.min(miny, data.coords[tets[i].v[j]][1]);
				minz = Math.min(minz, data.coords[tets[i].v[j]][2]);
				
				maxx = Math.max(maxx, data.coords[tets[i].v[j]][0]);
				maxy = Math.max(maxy, data.coords[tets[i].v[j]][1]);
				maxz = Math.max(maxz, data.coords[tets[i].v[j]][2]);
			}
			getArcs(tets[i]);
			for(int x = (int) minx;x < maxx; x ++) {
				for(int y = (int) miny;y < maxy;y ++) {
					for(int z = (int) minz;z < maxz;z ++) {
						double [] b = isInside(tets[i],x,y,z);
						if(b != null) {
							double dd = 0;
							for(int j = 0;j < 4;j ++) {
								dd += b[j] * data.fnVertices[tets[i].v[j]];
							}
							cl[x][y][z] = getColor(dd);
							ct[x][y][z] ++;
							gData[x][y][z] = dd;
						}
					}
				}
			}
		}
		
		for(int x = 0;x < size; x ++) {
			for(int y = 0;y < size;y ++) {
				for(int z = 0;z < size;z ++) {
					gData[x][y][z] += cl[x][y][z];
					max = Math.max(max, gData[x][y][z]);
				}
			}
		}
		System.out.println("Max Value : " + max);
	}

	private int getColor(double dd) {
		int ct = 0;
		for(Iterator<Integer> it = colV.iterator();it.hasNext();) {
			int v = it.next();
			if(data.fnVertices[v] >= dd) {
				return col.get(ct);
			}
			ct ++;
		}
		return col.get(ct - 1);
	}

	ArrayList<Integer> col = new ArrayList<Integer>(4);
	ArrayList<Integer> colarcs = new ArrayList<Integer>(4);
	ArrayList<Integer> colV = new ArrayList<Integer>(4);
	HashMap<Integer, Integer> prev = new HashMap<Integer, Integer>();
	ArrayList<Integer> q = new ArrayList<Integer>();
	HashSet<Integer> qq = new HashSet<Integer>();
	
	private void getArcs(Tetrahedron tet) {
		colarcs.clear();
		for(int i = 0;i < 4;i ++) {
			if(!colarcs.contains(arcs[tet.v[i]])) {
				colarcs.add(arcs[tet.v[i]]);
			}
		}
		int from = -1;
		int to = -1;
		for(Iterator<Integer> it = colarcs.iterator();it.hasNext();) {
			int a = it.next();
			if(from == -1) {
				from = rgData.arcs[a].from;
				to = rgData.arcs[a].to;
			} else {
				if(from > rgData.arcs[a].from) {
					from = rgData.arcs[a].from;
				}
				if(to < rgData.arcs[a].to) {
					to = rgData.arcs[a].to;
				}
			}
		}
		int v = from;
		col.clear();
		prev.clear();
		q.clear();
		qq.clear();
		qq.add(v);
		q.add(v);
		colV.clear();
		while(v != to && q.size() > 0) {
			v = q.remove(0);
			int barc = -1;
			for(Iterator<Integer> it = rgData.nodes[v].next.iterator();it.hasNext();) {
				int a = it.next();
				if(colarcs.contains(a)) {
					barc = a;
					break;
				}
			}			
			if(barc != -1) {
				int a = barc;
				Arc arc = rgData.arcs[a];
				if(arc.to <= to && !qq.contains(arc.to)) {
					q.add(arc.to);
					prev.put(arc.to, a);
				}
			} else {
				for(Iterator<Integer> it = rgData.nodes[v].next.iterator();it.hasNext();) {
					int a = it.next();
					Arc arc = rgData.arcs[a];
					if(arc.to <= to && !qq.contains(arc.to)) {
						q.add(arc.to);
						prev.put(arc.to, a);
					}
				}
			}
		}
		if(v != to) {
//			for(Iterator<Integer> it = colarcs.iterator();it.hasNext();) {
//				int a = it.next();
//				if(!colV.contains(rgData.arcs[a].from)) {
//					colV.add(rgData.arcs[a].from);
//				}
//				if(!colV.contains(rgData.arcs[a].to)) {
//					colV.add(rgData.arcs[a].to);
//				}
//			}
//			Collections.sort(colV);
			col.addAll(colarcs);
		} else {
			while(v != from) {
				int a = prev.get(v);
				col.add(0,a);
				v = rgData.arcs[a].from;
			}
		}
		
		for(Iterator<Integer> it = col.iterator();it.hasNext();) {
			int a = it.next();
			int vv = rgData.nodes[rgData.arcs[a].to].v;
			colV.add(vv);
		}
	}

	int noBytes = 1;
	double [] scale;
	int [] dim;
	double min, max;
	
	int [] cols;
	
	double [][] mat = new double[4][4];
	private double [] isInside(Tetrahedron t, int x, int y, int z) {
		mat[0][0] = data.coords[t.v[0]][0]; 
		mat[0][1] = data.coords[t.v[0]][1];
		mat[0][2] = data.coords[t.v[0]][2];
		mat[0][3] = 1;
		
		mat[1][0] = data.coords[t.v[1]][0]; 
		mat[1][1] = data.coords[t.v[1]][1];
		mat[1][2] = data.coords[t.v[1]][2];
		mat[1][3] = 1;
		
		mat[2][0] = data.coords[t.v[2]][0]; 
		mat[2][1] = data.coords[t.v[2]][1];
		mat[2][2] = data.coords[t.v[2]][2];
		mat[2][3] = 1;
		
		mat[3][0] = data.coords[t.v[3]][0]; 
		mat[3][1] = data.coords[t.v[3]][1];
		mat[3][2] = data.coords[t.v[3]][2];
		mat[3][3] = 1;

		double d0 = determinant(mat);
		
		mat[0][0] = x; 
		mat[0][1] = y;
		mat[0][2] = z;
		mat[0][3] = 1;
		double d1 = determinant(mat);
		
		
		mat[0][0] = data.coords[t.v[0]][0]; 
		mat[0][1] = data.coords[t.v[0]][1];
		mat[0][2] = data.coords[t.v[0]][2];
		mat[0][3] = 1;
		
		mat[1][0] = x; 
		mat[1][1] = y;
		mat[1][2] = z;
		mat[1][3] = 1;
		double d2 = determinant(mat);
		
		mat[1][0] = data.coords[t.v[1]][0]; 
		mat[1][1] = data.coords[t.v[1]][1];
		mat[1][2] = data.coords[t.v[1]][2];
		mat[1][3] = 1;
		
		mat[2][0] = x; 
		mat[2][1] = y;
		mat[2][2] = z;
		mat[2][3] = 1;
		double d3 = determinant(mat);
		
		mat[2][0] = data.coords[t.v[2]][0]; 
		mat[2][1] = data.coords[t.v[2]][1];
		mat[2][2] = data.coords[t.v[2]][2];
		mat[2][3] = 1;
		
		mat[3][0] = x; 
		mat[3][1] = y;
		mat[3][2] = z;
		mat[3][3] = 1;
		double d4 = determinant(mat);
		
		double d = d1 + d2 + d3 + d4;
		if(Math.abs(d - d0) >= 0.001) {
			System.out.println("Some Error");
			System.exit(0);
		}
		
		if(d == 0) {
			return null;
		}
		
		boolean inside = false;
		if(d0 >= 0 && d1 >= 0 && d2 >= 0 && d3 >= 0 && d4 >= 0) {
			inside = true;
		}
		if(d0 <= 0 && d1 <= 0 && d2 <= 0 && d3 <= 0 && d4 <= 0) {
			inside = true;
		}
		if(!inside) {
			return null;
		}
		double [] b = new double[4];
		b[0] = (double) (d1 / d);
		b[1] = (double) (d2 / d);
		b[2] = (double) (d3 / d);
		b[3] = (double) (d4 / d);
		return b;
	}
	
	public double determinant(double[][] mat) {

		double result = 0;

		if (mat.length == 1) {
			result = mat[0][0];
			return result;
		}

		if (mat.length == 2) {
			result = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
			return result;
		}

		for (int i = 0; i < mat[0].length; i++) {
			double temp[][] = new double[mat.length - 1][mat[0].length - 1];

			for (int j = 1; j < mat.length; j++) {
				for (int k = 0; k < mat[0].length; k++) {

					if (k < i) {
						temp[j - 1][k] = mat[j][k];
					} else if (k > i) {
						temp[j - 1][k - 1] = mat[j][k];
					}

				}
			}

			result += mat[0][i] * Math.pow(-1, (double) i) * determinant(temp);
		}
		return result;
	}


	boolean isGreater(int v1, int v2) {
		if (data.fnVertices[v1] > data.fnVertices[v2] || data.fnVertices[v1] == data.fnVertices[v2] && v1 > v2) {
			return true;
		}
		return false;
	}
	
	public static void main(String [] args) {
		
	}

	boolean raw = false;
	
	public void setRaw(boolean raw) {
		this.raw = raw;
	}

	void writeVRFFile(String op) {
		int maxx = size;
		int minx = 0;
		int maxy = size;
		int miny = 0;
		int maxz = size;
		int minz = 0;
		double scalex = 1;
		double scaley = 1;
		double scalez = 1;
		if(raw) {
			maxx = dim[0];
			maxy = dim[1];
			maxz = dim[2];
			
			scalex = scale[0];
			scaley = scale[1];
			scalez = scale[2];
		}
		
		int sizex = maxx - minx;
		int sizey = maxy - miny;
		int sizez = maxz - minz;
		System.out.println(sizex + " " + sizey + " " + sizez);
		try {
			RandomAccessFile outFile = new RandomAccessFile(op, "rw");
			outFile.writeShort(sizex);
			outFile.writeShort(sizey);
			outFile.writeShort(sizez);
			
			outFile.writeDouble(scalex);
			outFile.writeDouble(scaley);
			outFile.writeDouble(scalez);
			
			for(int z = minz;z < maxz;z ++) {
				for(int y = miny;y < maxy;y ++) {
					for(int x = minx;x < maxx;x ++) {
						outFile.writeDouble(gData[x][y][z]);
					}
				}
			}
			outFile.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void writeOffFile(String orig, String nu, int [] cols) {
		try {
			PrintStream p = new PrintStream(nu);
			FileInputStream ff = new FileInputStream(new File(orig));
			InputStreamReader reader = new InputStreamReader(ff);
			BufferedReader f = new BufferedReader(reader);
			String line = f.readLine();
			if(line.equalsIgnoreCase("OFF")) {
				line = f.readLine();
				p.println("OFF");
			}
			String [] s = Utilities.splitString(line);
			int nv = Integer.parseInt(s[0]);
			int nt = Integer.parseInt(s[1]);
			p.println(nv + " " + nt + " 0");
			for(int i = 0;i < nv;i ++) {
				line = f.readLine();
				p.println(line.trim() + " " + cols[i]);
			}
			for(int i = 0;i < nt;i ++) {
				line = f.readLine();
				p.println(line);
			}
			f.close();
			p.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void writePartition(String rawFile,int [] nodeComp, String opFile, boolean ascii) {
		try {
			if(ascii) {
				if(rawFile.endsWith(".dx")) {
					readDXFile(rawFile);
				} else {
					readTextFile(rawFile);
				}
				arcs = nodeComp;
				normalizeAndAssignArcs();
			} else {
				dim = new int[3];
				scale = new double[3];
				double [] data = GridReader.readVRF(rawFile, dim, scale);
				gData = new double[dim[0]][dim[1]][dim[2]];
				int ct = 0;
				arcs = nodeComp;
				min = Double.MAX_VALUE;
				max = -Double.MAX_VALUE;
				arcs = nodeComp;
				for(int k = 0;k < dim[2];k ++) {
					for(int j = 0;j < dim[1];j ++) {
						for(int i = 0;i < dim[0];i ++) {
							gData[i][j][k] = data[ct];
							min = Math.min(min, gData[i][j][k]);
							max = Math.max(max, gData[i][j][k]);
							gData[i][j][k] += arcs[ct ++];
						}
					}
				}
				System.out.println("MIN :: MAX ::::: " + min + " :: " + max);
			}
		} catch (Exception e) {
			e.printStackTrace();
			er("Error");
		}
		writeVRFFile(opFile);
	}

	private void normalizeAndAssignArcs() {
		int ct = 0;
		double diff = max - min;
		if(diff == 0) {
			diff = 1;
		}
		for(int k = 0;k < dim[2];k ++) {
			for(int j = 0;j < dim[1];j ++) {
				for(int i = 0;i < dim[0];i ++) {
					gData[i][j][k] = (gData[i][j][k] - min)/diff;
					if(gData[i][j][k] == 1) {
						gData[i][j][k] = 0.9999f;
					}
					gData[i][j][k] += arcs[ct ++];
				}
			}
		}
	}

	private int[] getDimension(String fileName) {
		try {
			int [] dim = new int[3];
			int in = fileName.lastIndexOf(".");
			fileName = fileName.substring(0,in);
			fileName += ".dim";
			FileInputStream ff = new FileInputStream(new File(fileName));
			InputStreamReader reader = new InputStreamReader(ff);
			BufferedReader f = new BufferedReader(reader);
			String s = f.readLine();
			String [] r = Utilities.splitString(s);
			dim[0] = Integer.parseInt(r[0].trim());
			dim[1] = Integer.parseInt(r[1].trim());
			dim[2] = Integer.parseInt(r[2].trim());
			f.close();
			return dim;
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	double [] getScale(String fileName) {
		try {
			double [] scale = new double [3];
			int in = fileName.lastIndexOf(".");
			fileName = fileName.substring(0,in);
			fileName += ".scale";
			FileInputStream ff = new FileInputStream(new File(fileName));
			InputStreamReader reader = new InputStreamReader(ff);
			BufferedReader f = new BufferedReader(reader);
			String s = f.readLine();
			String [] r = Utilities.splitString(s);
			scale[0] = Double.parseDouble(r[0].trim());
			scale[1] = Double.parseDouble(r[1].trim());
			scale[2] = Double.parseDouble(r[2].trim());
			f.close();
			return scale;
		} catch(Exception e) {
			double [] scale = new double[] {1,1,1};
			return scale;
		}
	}
	
	
	public void readRawFile(String fileName) throws Exception {
		dim = getDimension(fileName);
		scale = getScale(fileName);
		int nv = dim[0] * dim[1] * dim[2];
		RandomAccessFile fl = new RandomAccessFile(fileName,"r");
		pr("No. of Vertices : " + nv);
		gData = new double[dim[0]][dim[1]][dim[2]];
		min = Double.MAX_VALUE;
		max = -Double.MAX_VALUE;
		
		if(noBytes == 1) {
			// read vertices
			for(int k = 0;k < dim[2];k ++) {
				for(int j = 0;j < dim[1];j ++) {
					for(int i = 0;i < dim[0];i ++) {
						double fn;
						
						fn = fl.readUnsignedByte();
						
						double val = fn;
						gData[i][j][k] = val;
						min = Math.min(min, gData[i][j][k]);
						max = Math.max(max, gData[i][j][k]);
					}
				}
			}
			fl.close();
		} else {
			// read vertices
			
			for(int k = 0;k < dim[2];k ++) {
				for(int j = 0;j < dim[1];j ++) {
					for(int i = 0;i < dim[0];i ++) {
						double fn;
						
						int high = fl.readUnsignedByte() & 0xff;
						int low = fl.readUnsignedByte() & 0xff;

						// combine into a signed short.
						int data =  high << 8 | low;
						
						fn = (double) data;

						gData[i][j][k] = fn;
						min = Math.min(min, gData[i][j][k]);
						max = Math.max(max, gData[i][j][k]);
					}
				}
			}
			fl.close();
		}
	}
	
	public void readTextFile(String fileName) {
		try {
			min = Double.MAX_VALUE;
			max = -Double.MAX_VALUE;

			dim = getDimension(fileName);
			scale = getScale(fileName);
			FileInputStream ff = new FileInputStream(new File(fileName));
			InputStreamReader reader = new InputStreamReader(ff);
			BufferedReader f = new BufferedReader(reader);
			int nv = dim[0] * dim[1] * dim[2];
			pr("No. of Vertices : " + nv);
			gData = new double[dim[0]][dim[1]][dim[2]];
			for(int k = 0;k < dim[2];k ++) {
				for(int j = 0;j < dim[1];j ++) {
					for(int i = 0;i < dim[0];i ++) {
						String [] s = Utilities.splitString(f.readLine().trim());
						int l = s.length;
						double v = Double.parseDouble(s[l - 1].trim());
						
						gData[i][j][k] = v; 
						min = Math.min(min, gData[i][j][k]);
						max = Math.max(max, gData[i][j][k]);
					}
				}
			}
			f.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void readDXFile(String fileName) {
		try {
			dim = getDimension(fileName);
			scale = getScale(fileName);
			FileInputStream ff = new FileInputStream(new File(fileName));
			InputStreamReader reader = new InputStreamReader(ff);
			BufferedReader f = new BufferedReader(reader);
			int nv = dim[0] * dim[1] * dim[2];
			pr("No. of Vertices : " + nv);
			gData = new double[dim[0]][dim[1]][dim[2]];
			
			int ct = 0;
			while(ct != nv){
				String [] s = splitString(f.readLine().trim());
				for(int xx = 0;xx < s.length;xx ++) {
					double v = Double.parseDouble(s[xx].trim());
					int i = (ct % (dim[0] * dim[1])) % dim[0];
					int j = (ct % (dim[0] * dim[1])) / dim[0];
					int k = ct / (dim[0] * dim[1]);
					if(j == 1) {
						j = 1;
					}
					if(k == 1) {
						k = 1;
					}
					gData[i][j][k] = v; 
					min = Math.min(min, gData[i][j][k]);
					max = Math.max(max, gData[i][j][k]);
					ct ++;
				}
			}
			f.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
