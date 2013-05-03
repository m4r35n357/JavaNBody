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
	
	double ts = 0.005;
	
	int np = 0;
	
	List<Particle> particles = new ArrayList<Particle>();
	
	public Symplectic (int np) {		
		this.np = np;
	}
	
	public Symplectic (double g, double ts, int np) {
		this.g = g;
		this.ts = ts;
		this.np = np;
	}
	
	double distance (double x1, double y1, double z1, double x2, double y2, double z2) {
		return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2) + Math.pow(z2 - z1, 2));
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

	public void sympEuler (Coordinates first, Coordinates second) {
		first.update(this, 1.0);
		second.update(this, 1.0);
	}

	public void stormerVerlet2 (Coordinates first, Coordinates second) {
		first.update(this, 0.5);
		second.update(this, 1.0);
		first.update(this, 0.5);
	}

	public void stormerVerlet4 (Coordinates first, Coordinates second) {
		double c1, c2, c3, c4, d1, d2, d3;
		double cubeRoot2 = Math.pow(2.0, 1.0 / 3.0);
		double denom = 2.0 - cubeRoot2;
		c1 = c4 = 1.0 / (2.0 * denom);
		c2 = c3 = (1.0 - cubeRoot2) / (2.0 * denom);
		d1 = d3 = 1.0 / denom;
		d2 = - cubeRoot2 / denom;
		first.update(this, c1);
		second.update(this, d1);
		first.update(this, c2);
		second.update(this, d2);
		first.update(this, c3);
		second.update(this, d3);
		first.update(this, c4);
	}

	public void setParticles (List<Particle> particles) {
		this.particles = particles;
	}

	/**
	 * @param args
	 * @throws InterruptedException 
	 */
	public static void main (String[] args) throws InterruptedException {
		boolean debug = true;
		boolean update;
		long n = 0;
		double H0, Hcurrent, Hmin, Hmax, error, snr;
		List<Particle> particles = new ArrayList<Particle>();
		// 2-body test
//		particles.add(new Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0));
//		particles.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		// 8-body test
		particles.add(new Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 100.0));
		particles.add(new Particle(0.0, 1.5, 0.4, -3.4, 0.0, -0.2, 2.0));
		particles.add(new Particle(-2.0, 0.0, -0.4, 0.0, -4.0, 0.2, 3.0));
		particles.add(new Particle(3.0, 0.0, -0.2, 0.0, 4.0, -0.1, 5.0));
		particles.add(new Particle(0.0, -4.0, 0.1, 4.6, 0.0, 0.1, 4.0));
		particles.add(new Particle(-4.0, 0.0, -0.1, 0.0, -2.8, -0.2, 3.0));
		particles.add(new Particle(2.0, 0.0, -0.3, 0.0, 4.4, 0.2, 3.0));
		particles.add(new Particle(0.0, 3.0, -0.2, -5.0, 0.0, -0.1, 4.0));
		Symplectic symplectic = new Symplectic(particles.size());
		symplectic.setParticles(particles);
		H0 = symplectic.hamiltonian();
		Hmin = H0;
		Hmax = H0;
		error = 0.0;
		while (n <= 10000) {
			symplectic.stormerVerlet4(Coordinates.POSITION, Coordinates.MOMENTUM);
			if (debug) {
				update = false;
				Hcurrent = symplectic.hamiltonian();
				if (Hcurrent < Hmin) {
					Hmin = Hcurrent;
					error += Math.abs(Hcurrent - H0);
					update = true;
				} else if (Hcurrent > Hmax) {
					Hmax = Hcurrent;
					error += Math.abs(Hcurrent - H0);
					update = true;
				}
//				if (update && ((n % 1000) == 0)) {
				if ((n % 1000) == 0) {
					snr = 10.0 * Math.log10(Math.abs(H0 / error));
					System.out.printf("n: %9d, Hamiltonian: %.9e, Start: %.9e, Hmin %.9e, Hmax %.9e, Error: %.9e, SNR: %6.1f%n", n, symplectic.hamiltonian(), H0, Hmin, Hmax, error, snr);
				}
			}
			n += 1;
		}
	}

}
