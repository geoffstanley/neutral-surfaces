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
package vgl.iisc.utils;

public class Pair<A, B> {
    private A first;
    private B second;

    public Pair(A first, B second) {
        super();
        this.first = first;
        this.second = second;
    }

    public int hashCode() {
        int hashFirst = first != null ? first.hashCode() : 0;
        int hashSecond = second != null ? second.hashCode() : 0;

        return (hashFirst + hashSecond) * hashSecond + hashFirst;
    }

    @SuppressWarnings("unchecked")
	public boolean equals(Object other) {
        if (other instanceof Pair) {
                Pair<A, B> otherPair = (Pair<A, B>) other;
                return 
                ((  this.first == otherPair.first ||
                        ( this.first != null && otherPair.first != null &&
                          this.first.equals(otherPair.first))) &&
                 (      this.second == otherPair.second ||
                        ( this.second != null && otherPair.second != null &&
                          this.second.equals(otherPair.second))) );
        }

        return false;
    }

    public String toString()
    { 
           return "(" + first.toString() + ", " + second.toString() + ")"; 
    }

    public A getFirst() {
        return first;
    }

    public void setFirst(A first) {
        this.first = first;
    }

    public B getSecond() {
        return second;
    }

    public void setSecond(B second) {
        this.second = second;
    }
}

