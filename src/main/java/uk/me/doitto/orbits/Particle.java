package uk.me.doitto.orbits;

import java.io.Serializable;

/**
 * @author ian
 * <p>
 * Holds per-particle data
 */
public class Particle implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	double qX, qY, qZ, pX, pY, pZ, mass;
	
	public Particle (double qX, double qY, double qZ, double pX, double pY, double pZ, double mass) {
		this.qX = qX;
		this.qY = qY;
		this.qZ = qZ;
		this.pX = pX;
		this.pY = pY;
		this.pZ = pZ;
		this.mass = mass;
	}
	
	/**
	 * Writes out current particle state as JSON-formatted text
	 * @return the JSON string
	 */
	public String toString () {
		return String.format("{\"qX\":%.18f,\"qY\":%.18f,\"qZ\":%.18f,\"pX\":%.18f,\"pY\":%.18f,\"pZ\":%.18f,\"mass\": %3.1f}", qX, qY, qZ, pX, pY, pZ, mass);
//		return "{\"qX\":" + qX + ",\"qY\":" + qY + ",\"qZ\":" + qZ + ",\"pX\":" + pX + ",\"pY\":" + pY + ",\"pZ\":" + pZ + ",\"mass\":" + mass + "}";
	}
}
