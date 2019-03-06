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
package vgl.iisc.external.loader;

import java.io.File;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * This class is used to track and select the appropriate loader to load the input mesh. 
 * 
 */
public class DataLoader {

	private static DataLoader loader;

	private HashMap<String, String> loaderMap;
	
	private DataLoader() {
		loaderMap = new HashMap<String, String>();
		
		try {
			File file = new File("loaders.xml");
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(file);
			doc.getDocumentElement().normalize();
			String root = doc.getDocumentElement().getNodeName();
			if (!root.trim().equals("loaders")) {
				System.err.println("Invalid loader file.");
				System.exit(0);
			}

			NodeList meshList = doc.getElementsByTagName("mesh");

			for (int s = 0; s < meshList.getLength(); s++) {
				Node node = meshList.item(s);
				if (node.getNodeType() == Node.ELEMENT_NODE) {
					Element mesh = (Element) node;
					NodeList typeList = mesh.getElementsByTagName("type");
					Element type = (Element) typeList.item(0);
					NodeList format = type.getChildNodes();
					String typeName = ((Node) format.item(0)).getNodeValue().trim();
					NodeList classList = mesh.getElementsByTagName("class");
					Element classElmnt = (Element) classList.item(0);
					NodeList classType = classElmnt.getChildNodes();
					String className = ((Node) classType.item(0)).getNodeValue().trim();
					
					loaderMap.put(typeName, className);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	private MeshLoader getMeshLoader(String type) {
		String className = loaderMap.get(type.trim());
		if(className == null) {
			System.err.println("Class for type " + type + " not defined in loaders.xml");
			System.exit(0);
		}
		try {
			Class<?> loaderClass = Class.forName(className);
			Object loaderObj = loaderClass.newInstance();
			if(loaderObj instanceof MeshLoader) {
				return (MeshLoader) loaderObj;
			} else {
				System.err.println("Class " + className + " does not extend class iisc.vgl.external.loader.MeshLoader");
				System.exit(0);
			}
		} catch (ClassNotFoundException e) {
			System.err.println("Class" + className + " for type " + type + " not found. Check if it is in added in the program's classpath.");
			System.exit(0);
		} catch (InstantiationException e) {
			System.err.println("Unable to instantiate Class" + className + " for type " + type + ".");
			e.printStackTrace();
			System.exit(0);
		} catch (IllegalAccessException e) {
			System.err.println("Unable to instantiate Class" + className + " for type " + type + ".");
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}
	
	/**
	 * This method is called by different classes in order to obtain the loader to load various file formats.
	 * 
	 * @param type Specifies the input format. Note that it should have an entry in the <i>loaders.xml</i> file. 
	 * @return The custom loader corresponding to the given type 
	 */
	public static MeshLoader getLoader(String type) {
		if (loader == null) {
			loader = new DataLoader();
		}
		return loader.getMeshLoader(type);
	}
}
