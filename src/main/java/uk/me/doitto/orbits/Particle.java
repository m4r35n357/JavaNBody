/**
 * 
 */
package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
public class Particle {
	
	public double qX, qY, qZ, pX, pY, pZ, mass;
	
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

}
