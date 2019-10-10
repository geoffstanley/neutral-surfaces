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
import java.util.StringTokenizer;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;
import vgl.iisc.external.types.Triangle;
import vgl.iisc.external.types.Vertex;

/**
 * 
 * This class is used to load a tetrahedral mesh. It assumes the following format for the mesh:
 * <ul>
 * <li> An optional first line with "OFF" specifying file type 
 * <li>	The next line specifies the no. of vertices (nv) followed by the number of triangles (nv) (space seperated) </li>
 * <li>	The next nv lines contains <br/>
 * 		x y z [f] <br/>
 * 		where x, y & z specify the co-ordinates of the vertex and f specifies the function value. </li>
 * <li>	the next nt lines has <br/>
 * 		[3] v1 v2 v3 <br/>
 * 		where v1, v2 and v3 are the vertex indices of the vertices that form the triangles (the 3 is optional) </li>
 * </ul>
 * 
 */
public class OffLoader implements MeshLoader {

	private BufferedReader reader;
	private int noVertices;
	private int noTris;
	private int curVertex;
	private int curTri;
	private String mesh;
	
	@Override
	public void setInputFile(String inputMesh) {
		try {
			mesh = inputMesh;
			reader = new BufferedReader(new FileReader(inputMesh));
			String s = reader.readLine();
			if(s.trim().equalsIgnoreCase("OFF")) {
				s = reader.readLine();
			}
			String[] r = splitString(s);
			noVertices = Integer.parseInt(r[0].trim());
			if (r.length == 1) {
				s = reader.readLine();
				r = splitString(s);
				noTris = Integer.parseInt(r[0].trim());
			} else {
				noTris = Integer.parseInt(r[1].trim());
			}

			curVertex = 0;
			curTri = 0;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	@Override
	public int getVertexCount() {
		return noVertices;
	}

	@Override
	public int getSimplexCount() {
		return noTris;
	}
	
	@Override
	public Simplex getNextSimplex() {
		try {
			if (curVertex < noVertices) {
				double x, y, z, fn;
				String s = reader.readLine();
				if(s.trim().equalsIgnoreCase("")) {
					s = reader.readLine();
				}
				String[] r = splitString(s);
				x = Double.parseDouble(r[0].trim());
				y = Double.parseDouble(r[1].trim());
				z = Double.parseDouble(r[2].trim());
				if (r.length == 4 && r[3].trim().length() > 0) {
					fn = Double.parseDouble(r[3].trim());
				} else {
					fn = -1;
				}
				Vertex vertex = new Vertex();
				vertex.c = new double[3];
				vertex.c[0] = x;
				vertex.c[1] = y;
				vertex.c[2] = z;
				vertex.f = fn;

				curVertex++;
				return vertex;
			}
			if (curTri < noTris) {
				String s = reader.readLine();
				String[] r = splitString(s);
				int v1 = -1;
				int v2 = -1;
				int v3 = -1;
				if(r.length >= 4) {
					int f = Integer.parseInt(r[0]);
					if(f != 3) {
						System.err.println("Invalid input : No. of vertices in a face is not 3");
						System.exit(0);
					}
					v1 = Integer.parseInt(r[1]);
					v2 = Integer.parseInt(r[2]);
					v3 = Integer.parseInt(r[3]);
				} else if(r.length == 3) {
					v1 = Integer.parseInt(r[0]);
					v2 = Integer.parseInt(r[1]);
					v3 = Integer.parseInt(r[2]);
				} else {
					System.err.println("Invalid input");
					System.exit(0);
				}

				Triangle tri = new Triangle();
				tri.v1 = v1;
				tri.v2 = v2;
				tri.v3 = v3;

				curTri++;
				return tri;
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
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
		return null;
	}
}
