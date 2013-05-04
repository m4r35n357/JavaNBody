/**
 * 
 */
package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
public interface IIntegrators {
	void solve (Symplectic s, Coordinates first, Coordinates second);
}
