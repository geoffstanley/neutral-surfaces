/*
 *	Copyright (C) 2010 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of libRG, a library to compute Reeb graphs.
 *
 *	libRG is free software: you can redistribute it and/or modify
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
 *	Modified by : -- 
 *	Date : --
 *	Changes  : --
 */
package vgl.iisc.external.loader.impl;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.StringTokenizer;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;
import vgl.iisc.external.types.Triangle;
import vgl.iisc.external.types.Vertex;


/**
 * 
 * This class is used to load an arbitrary simplicial complex. It assumes the following format for the mesh:
 * <ul>
 * <li>	First line specifies the dimension (d) of the input </li>
 * <li>	The next line specifies the no. of vertices (nv) followed by the number of simplices (ns) (space seperated) </li>
 * <li>	The next nv lines contains <br/>
 * 		c<sub>1</sub> c<sub>2</sub> ... c<sub>d</sub> [f] <br/>
 * 		where c<sub>i</sub> specifies the i<sup>th</sup> co-ordinate of the vertex and f specifies the function value. </li>
 * <li>	the next ns lines has <br/>
 * 		(l + 1) v<sub>1</sub> v<sub>2</sub> ... v<sub>l+1</sub> <br/>
 * 		where l is the dimension of the simplex and v<sub>i</sub> is the index of the i<sup>th</sup> vertex of the simplex. </li>
 * </ul>
 * 
 */
public class SimLoader implements MeshLoader {
	
	class TmpVertex {
		Vertex v;
		int index;
	}
	
	class TmpSimplex {
		int [] v;
	}
	
	private BufferedReader reader;
	private int noVertices;
	private int noSims;
	private int curVertex;
	private int curSim;
	private int dim;
	private String mesh;

	private ArrayList<Triangle> tris = new ArrayList<Triangle>();
	private ArrayList<TmpVertex> verts = new ArrayList<TmpVertex>();
	private ArrayList<TmpVertex> sortVerts = new ArrayList<TmpVertex>();
	private ArrayList<TmpSimplex> sims = new ArrayList<TmpSimplex>();
	
	private int [] vMap;
	
	@Override
	public void setInputFile(String inputMesh) {
		try {
			mesh = inputMesh;
			reader = new BufferedReader(new FileReader(inputMesh));
			String s = reader.readLine();
			dim = Integer.parseInt(s.trim());
			s = reader.readLine();
			String[] r = splitString(s);
			noVertices = Integer.parseInt(r[0].trim());
			noSims = Integer.parseInt(r[1].trim());

			System.out.println("No. of Vertices : " + noVertices);
			System.out.println("No. of Simplices : " + noSims);

			curVertex = 0;
			curSim = 0;
			
			for(int n = 0;n < noVertices;n ++) {
				double fn;
				s = reader.readLine();
				r = splitString(s);
				if(r.length < dim) {
					System.out.println("Invalid d-dimensional point");
					System.exit(0);
				}
				Vertex vertex = new Vertex();
				vertex.c = new double[dim];
				for(int i = 0;i < dim;i ++) {
					vertex.c[i] = Double.parseDouble(r[i].trim());
				}
				if (r.length == dim+1 && r[dim].trim().length() > 0) {
					fn = Double.parseDouble(r[dim].trim());
				} else {
					fn = -1;
				}
				vertex.f = fn;

				TmpVertex tv = new TmpVertex();
				tv.v = vertex;
				tv.index = n;
				verts.add(tv);
			}
			
			for(int n = 0;n < noSims;n ++) {
				s = reader.readLine();
				r = splitString(s);
				int l = Integer.parseInt(r[0]); 
				int [] v = new int [l];
				
				for(int i = 0;i < l;i ++) {
					v[i] = Integer.parseInt(r[i+1].trim());
				}
				TmpSimplex sim = new TmpSimplex();
				if(v.length > 2) {
					sim.v = v;
				} else if (v.length == 2) {
					int [] vv = new int[3];
					TmpVertex tv = new TmpVertex();
					TmpVertex tv1 = verts.get(v[0]);
					TmpVertex tv2 = verts.get(v[1]);
					if(tv2.index > tv2.index) {
						TmpVertex tmp = tv2;
						tv2 = tv1;
						tv1 = tmp;
					}
					tv.v = new Vertex();
					tv.v.c = new double[tv1.v.c.length];
					for(int d = 0;d < tv1.v.c.length; d ++) {
						tv.v.c[d] = (tv1.v.c[d] + tv2.v.c[d]) / 2;
					}
					tv.v.f = (tv1.v.f + tv2.v.f) / 2;
					tv.index = tv2.index;
					tv2.index = verts.size();
					verts.add(tv);
					vv[0] = v[0];
					vv[1] = v[1];
					vv[2] = verts.size() - 1;
					sim.v = vv;
				} else {
					System.out.println("Invalid simplex");
					System.exit(0);
				}
				sims.add(sim);
			}
			sortVerts.addAll(verts);
			Collections.sort(sortVerts, new Comparator<TmpVertex>() {

				@Override
				public int compare(TmpVertex o1, TmpVertex o2) {
					return o1.index - o2.index;
				}
			});
			vMap = new int[sortVerts.size()];
			Arrays.fill(vMap, -1);
			for(int i = 0;i < noVertices;i ++) {
				vMap[verts.get(i).index] = i;
			}
			noVertices = sortVerts.size();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	@Override
	public int getVertexCount() {
		return verts.size();
	}

	@Override
	public int getSimplexCount() {
		return noSims;
	}

	@Override
	public Simplex getNextSimplex() {
		try {
			if (curVertex < noVertices) {
				return sortVerts.get(curVertex ++).v;
			}
			if(!tris.isEmpty()) {
				 return tris.remove(0);
			}
			if (curSim < noSims) {
				TmpSimplex sim = sims.get(curSim);
				int l = sim.v.length;
				int [] v = new int [l];
				
				for(int i = 0;i < l;i ++) {
					v[i] = verts.get(sim.v[i]).index;
				}
				tris.clear();
				generateTris(v, -1, -1);
				curSim++;
				return tris.remove(0);
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}

	private void generateTris(int[] v, int v1, int v2) {
		if(v1 == -1) {
			for(int i = 0;i < v.length - 2;i ++) {
				generateTris(v, i, v2);
			}
		} else if(v2 == -1) {
			for(int i = v1 + 1;i < v.length - 1;i ++) {
				generateTris(v, v1, i);
			}
		} else {
			for(int i = v2 + 1;i < v.length;i ++) {
				Triangle t = new Triangle();
				t.v1 = v[v1];
				t.v2 = v[v2];
				t.v3 = v[i];
				tris.add(t);
			}
		}
	}

	private static String[] splitString(String s) {
		String[] ret = null;
		StringTokenizer tok = new StringTokenizer(s);
		ret = new String[tok.countTokens()];
		for (int i = 0; i < ret.length; i++) {
			ret[i] = tok.nextToken();
		}
		return ret;
	}
	
	@Override
	public void reset() {
		try {
			reader.close();
			setInputFile(mesh);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	@Override
	public int[] getVertexMap() {
		return vMap;
	}

}
