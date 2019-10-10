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


import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.StringTokenizer;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;
import vgl.iisc.external.types.Triangle;
import vgl.iisc.external.types.Vertex;

/**
 * 
 * This class is used to load a .ply surface mesh. It assumes the following format for the mesh:
 * 
 */
@SuppressWarnings("deprecation")
public class PlyLoader implements MeshLoader {
	
	public static enum DataType {CHAR, UCHAR, SHORT, USHORT, INT, UINT, FLOAT, DOUBLE};
	
//	private BufferedReader reader;
	DataInputStream bufReader;
	private int noVertices;
	private int noTris;
	private int curVertex;
	private int curTri;
	private String mesh;
	boolean ascii;
	boolean bigEndian;
	DataType xType, yType, zType;
	byte [] buf;
	int noProps;
	
	@Override
	public void setInputFile(String inputMesh) {
		try {
			mesh = inputMesh;
			bufReader = new DataInputStream(new BufferedInputStream(new FileInputStream(inputMesh)));
			String s = bufReader.readLine();
			if(!s.trim().equalsIgnoreCase("ply")) {
				System.err.println("Not a ply file");
				System.exit(0);
			}
			boolean vertexProp = false;
			int no = 0;
			while(!s.equalsIgnoreCase("end_header")) {
				System.out.println(s);
				s = bufReader.readLine();
				String[] r = splitString(s);
				if(r.length > 1 && r[0].trim().equalsIgnoreCase("element")){
					if(r[1].trim().equalsIgnoreCase("vertex")) {
						noVertices = Integer.parseInt(r[2].trim());
						vertexProp = true;
					}
					if(r[1].trim().equalsIgnoreCase("face")) {
						noTris = Integer.parseInt(r[2].trim());
						vertexProp = false;
					}
				}
				if(r.length > 1 && r[0].trim().equalsIgnoreCase("format")){
					if(r[1].trim().equalsIgnoreCase("ascii")) {
						System.out.println("Ascii format");
						ascii = true;
					} else {
						System.out.println("Binary format : ");
						if(r[1].trim().equals("binary_big_endian")) {
							System.out.println("Big endian");
							bigEndian = true;
						} else {
							System.out.println("little endian");
							bigEndian = false;
							System.err.println("little endian currently not supported");
							System.exit(0);
						}
						ascii = false;
					}
					System.out.println("Version : " + r[2].trim());
				}
				if(r.length > 1 && r[0].trim().equalsIgnoreCase("property")){
					if(vertexProp) {
						no ++;
					}
					if(r[2].trim().equalsIgnoreCase("x")) {
						if(r[1].equalsIgnoreCase("double")) {
							xType = DataType.FLOAT;
						} else {
							System.err.println("Data type currently not supported");
							System.exit(0);
						}
					}
					if(r[2].trim().equalsIgnoreCase("y")) {
						if(r[1].equalsIgnoreCase("double")) {
							yType = DataType.FLOAT;
						} else {
							System.err.println("Data type currently not supported");
							System.exit(0);
						}
					}
					if(r[2].trim().equalsIgnoreCase("z")) {
						if(r[1].equalsIgnoreCase("double")) {
							zType = DataType.FLOAT;
						} else {
							System.err.println("Data type currently not supported");
							System.exit(0);
						}
					}
					if(r[r.length - 1].trim().equalsIgnoreCase("vertex_indices")) {
						System.out.println(s);
					}
				}
			}
			
			System.out.println("No. of Vertices : " + noVertices);
			System.out.println("No. of Triangles : " + noTris);
//			System.out.println("No. of properties per vertex : " + no);
			noProps = no;
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
				if(ascii) {
					String s = bufReader.readLine();
					if(s.trim().equalsIgnoreCase("")) {
						s = bufReader.readLine();
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
				} else {
					x = readDouble();
					y = readDouble();
					z = readDouble();
					if(noProps > 3) {
						for(int _x = 3;_x < noProps;_x ++) {
							readDouble();	
						}
					}

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
				int v1 = -1;
				int v2 = -1;
				int v3 = -1;
				if(ascii) {
					String s = bufReader.readLine();
					String[] r = splitString(s);
					if(r.length == 4) {
						int f = Integer.parseInt(r[0]);
						if(f != 3) {
							System.err.println("Invalid input : No. of vertices in a face is not 3");
							System.exit(0);
						}
						v1 = Integer.parseInt(r[1]);
						v2 = Integer.parseInt(r[2]);
						v3 = Integer.parseInt(r[3]);
//					} else if(r.length == 3) {
//						v1 = Integer.parseInt(r[0]);
//						v2 = Integer.parseInt(r[1]);
//						v3 = Integer.parseInt(r[2]);
					} else {
						System.err.println("Invalid input");
						System.exit(0);
					}
				} else {
					int nf = readByte();
					if(nf != 3) {
						System.err.println("Only triangular faces are currently supported");
						System.exit(0);
					}
					v1 = readInt();
					v2 = readInt();
					v3 = readInt();
				}

				Triangle tri = new Triangle();
				tri.v1 = v1;
				tri.v2 = v2;
				tri.v3 = v3;

				curTri++;
				return tri;
			}
			bufReader.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}

	private double readDouble() throws IOException {
		return bufReader.readDouble();
	}

	private int readInt() throws IOException {
		return bufReader.readInt();
	}
	
	private int readByte() throws IOException {
		return bufReader.readUnsignedByte();
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
			bufReader.close();
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
