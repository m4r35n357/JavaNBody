/**
 * 
 */
package uk.me.doitto.orbits;

import java.io.Serializable;

/**
 * @author ian
 *
 */
public class Particle implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private static double DENSITY = 1.0e9;
	
	double qX, qY, qZ, pX, pY, pZ, mass;
	
	/**
	 * 
	 */
	public Particle (double qX, double qY, double qZ, double pX ,double pY ,double pZ, double mass) {
		this.qX = qX;
		this.qY = qY;
		this.qZ = qZ;
		this.pX = pX;
		this.pY = pY;
		this.pZ = pZ;
		this.mass = mass;
	}

	public double radius () {
		return Math.pow(mass / DENSITY, 1.0 / 3.0);
	}
	
}
