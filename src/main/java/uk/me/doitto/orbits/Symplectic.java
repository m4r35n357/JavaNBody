/**
 * 
 */
package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.List;

/**
 * @author ian
 *
 */
public class Symplectic {
	
	double g = 0.05;
	
	double ts = 0.001;
	
	int np = 0;
	
	List<Particle> particles = new ArrayList<Particle>();
	
	public Symplectic (InitialConditions ic) {
		this.particles = ic.bodies;
		this.np = ic.bodies.size();
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
		boolean debug = true;
		long n = 0;
		double H0, Hcurrent, Hmin, Hmax, error, snr;
		Symplectic s = new Symplectic(InitialConditions.EIGHT_BODY);
		H0 = s.hamiltonian();
		Hmin = H0;
		Hmax = H0;
		error = 0.0;
		while (n <= 400000) {
			Integrators.STORMER_VERLET_4.solve(s, Coordinates.Q, Coordinates.P);
			if (debug) {
				Hcurrent = s.hamiltonian();
				if (Hcurrent < Hmin) {
					Hmin = Hcurrent;
					error += Math.abs(Hcurrent - H0);
				} else if (Hcurrent > Hmax) {
					Hmax = Hcurrent;
					error += Math.abs(Hcurrent - H0);
				}
				if ((n % 1000) == 0) {
					snr = 10.0 * Math.log10(Math.abs((H0 - Hcurrent) / H0));
					System.out.printf("n: %9d, Hamiltonian: %.9e, Start: %.9e, Hmin %.9e, Hmax %.9e, Error: %.9e, ER: %6.1f%n", n, s.hamiltonian(), H0, Hmin, Hmax, error, snr);
				}
			}
			n += 1;
		}
	}
}
