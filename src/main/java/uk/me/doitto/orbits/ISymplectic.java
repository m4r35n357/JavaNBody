/**
 * 
 */
package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
interface ISymplectic {
	
	/**
	 * Total (kinetic + potential) energy of the system
	 * @return the total energy
	 */
	public double hamiltonian ();
	
	/**
	 * Position update implements dH/dp, which in this case is a function of p only
	 * @param c composition coefficient
	 */
	void updateQ (double c);
	
	/**
	 * Momentum update implements -dH/dq, which in this case is a function of q only
	 * @param c composition coefficient
	 */
	void updateP (double c);
}
