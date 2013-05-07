/**
 * 
 */
package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.List;
import static uk.me.doitto.orbits.Coordinates.*;
import static uk.me.doitto.orbits.Integrators.*;
import static uk.me.doitto.orbits.InitialConditions.*;

/**
 * @author ian
 *
 */
public class Symplectic {
	
	double g = 0.0;
	
	double ts = 0.0;
	
	int np = 0;
	
	List<Particle> particles = new ArrayList<Particle>();
	
	public Symplectic (double g, double ts, InitialConditions ic) {
		this.particles = ic.bodies;
		this.np = ic.bodies.size();
		this.g = g;
		this.ts = ts;
	}
	
	public Symplectic(double g, double ts, List<Particle> bodies) {
		// destroy reference to input array so the client can't change i
		this.particles = new ArrayList<Particle>(bodies);
		this.np = bodies.size();
		this.g = g;
		this.ts = ts;
	}

	public List<Particle> getParticles() {
		// return a copy
		return new ArrayList<Particle>(particles);
	}

	double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
		return Math.sqrt(Math.pow(xB - xA, 2) + Math.pow(yB - yA, 2) + Math.pow(zB - zA, 2));
	}
	
	public double hamiltonian () {
		Particle a, b;
		int i, j;
		double totalEnergy = 0.0;
		for (i = 0; i < np; i++) {
			a = particles.get(i);
			totalEnergy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass;
			for (j = 0; j < np; j++) {
				if (i > j) {
					b = particles.get(j);
					totalEnergy -= g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
				}
			}
		}
		return totalEnergy;
	}

	/**
	 * @param args
	 */
	public static void main (String[] args) {
		double h0, hMin, hMax;
		boolean debug = true;
		long n = 0;
		Symplectic s = new Symplectic(0.05, 0.001, EIGHT_BODY);
		h0 = s.hamiltonian();
		hMin = h0;
		hMax = h0;
		while (n <= 400000) {
			STORMER_VERLET_4.solve(s, Q, P);
			if (debug) {
				double hNow = s.hamiltonian();
				double dH = hNow - h0;
				if (hNow < hMin) {
					hMin = hNow;
				} else if (hNow > hMax) {
					hMax = hNow;
				}
				if ((n % 1000) == 0) {
					System.out.printf("n: %9d, Hamiltonian: %.9e, Start: %.9e, Hmin %.9e, Hmax %.9e, Error: %.3e, ER: %6.1f%n", n, hNow, h0, hMin, hMax, Math.abs(dH), 10.0 * Math.log10(Math.abs(dH / h0)));
				}
			}
			n += 1;
		}
	}
}
