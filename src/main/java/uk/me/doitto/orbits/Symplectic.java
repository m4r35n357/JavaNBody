package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.List;
import static uk.me.doitto.orbits.PhaseSpace.*;
import static uk.me.doitto.orbits.Integrator.*;
import static uk.me.doitto.orbits.Scenario.*;

/**
 * @author ian
 *
 * Main class for symplectic integrator simulations
 */
public class Symplectic {
	
	double iterations = 0.0;
	
	double g = 0.0;
	
	double timeStep = 0.0;
	
	int np = 0;
	
	List<Particle> particles = new ArrayList<Particle>();
	
	Integrator integrator;
	
	public Symplectic (Scenario ic, Integrator integrator) {
		this.particles = ic.bodies;
		this.np = ic.bodies.size();
		this.g = ic.g;
		this.timeStep = ic.ts;
		this.iterations = ic.simulationTime / ic.ts;
		this.integrator = integrator;
	}
	
	public Symplectic (double g, double simulationTime, double timeStep, List<Particle> bodies, String integrator) {
		// destroy reference to input array so the client can't change i
		this.particles = new ArrayList<Particle>(bodies);
		this.np = bodies.size();
		this.g = g;
		this.timeStep = timeStep;
		this.iterations = simulationTime / timeStep;
		this.integrator = Integrator.valueOf(integrator);
	}

	public List<Particle> getParticles () {
		// return a copy
		return new ArrayList<Particle>(particles);
	}

	double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
		return Math.sqrt(Math.pow(xB - xA, 2) + Math.pow(yB - yA, 2) + Math.pow(zB - zA, 2));
	}
	
	public double hamiltonian () {
		Particle a, b;
		double totalEnergy = 0.0;
		for (int i = 0; i < np; i++) {
			a = particles.get(i);
			totalEnergy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass;
			for (int j = 0; j < np; j++) {
				if (i > j) {
					b = particles.get(j);
					totalEnergy -= g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
				}
			}
		}
		return totalEnergy;
	}

	/**
	 * Test method for symplectic integrators
	 * 
	 * @param args None defined
	 */
	public static void main (String[] args) {
		Symplectic s = new Symplectic(THREE_BODY, STORMER_VERLET_4);
		long n = 0;
		double h0 = s.hamiltonian();
		double hMin = h0;
		double hMax = h0;
		while (n <= s.iterations) {
			s.integrator.solve(s, Q, P);
			double hNow = s.hamiltonian();
			double dH = hNow - h0;
			if (hNow < hMin) {
				hMin = hNow;
			} else if (hNow > hMax) {
				hMax = hNow;
			}
			if ((n % 1000) == 0) {
				StringBuilder json = new StringBuilder("[");
				for (Particle p : s.particles) {
					json.append("{\"Qx\":" + p.qX + ",\"Qy\":" + p.qY + ",\"Qz\":" + p.qZ + ",\"Px\":" + p.pX + ",\"Py\":" + p.pY + ",\"Pz\":" + p.pZ + "},");
				}
				System.out.println(json + "]");
				System.out.printf("t:%7.0f, H: %.9e, H0: %.9e, H-: %.9e, H+: %.9e, E: %.1e, ER: %6.1f dBH%n", n * s.timeStep, hNow, h0, hMin, hMax, Math.abs(dH), 10.0 * Math.log10(Math.abs(dH / h0)));
			}
			n += 1;
		}
	}
}
