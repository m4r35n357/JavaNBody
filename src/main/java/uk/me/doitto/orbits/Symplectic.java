package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author ian
 * <p>
 * Top level class for symplectic integrator simulations
 */
public class Symplectic {
	
	private final static double OUTPUT_TIME_GRANULARITY = 0.01;
	
	public final double iterations, g, timeStep, errorLimit;
	
	public final int outputInterval, np;
	
	public final List<Particle> particles;
	
	private final Integrator integrator;
	
	/**
	 * For calling from Java code
	 */
	public Symplectic (Scenario ic, Integrator integrator) {
		this.particles = ic.bodies;
		this.np = ic.bodies.size();
		this.g = ic.g;
		this.timeStep = ic.ts;
		this.errorLimit = ic.errorLimit;
		this.outputInterval = (int) Math.round(OUTPUT_TIME_GRANULARITY / ic.ts);
		this.iterations = ic.simulationTime / ic.ts;
		this.integrator = integrator;
	}
	
	/**
	 * For calling from other languages, or remotes
	 */
	public Symplectic (double g, double simulationTime, double timeStep, double errorLimit, List<Particle> bodies, String integrator) {
		this.particles = new ArrayList<Particle>(bodies);
		this.np = bodies.size();
		this.g = g;
		this.timeStep = timeStep;
		this.errorLimit = errorLimit;
		this.outputInterval = (int) Math.round(OUTPUT_TIME_GRANULARITY / timeStep);
		this.iterations = simulationTime / timeStep;
		this.integrator = Integrator.valueOf(integrator);
	}

	/**
	 * Separation between points A and B
	 * @param xA x coordinate of point A
	 * @param yA y coordinate of point A
	 * @param zA z coordinate of point A
	 * @param xB x coordinate of point B
	 * @param yB y coordinate of point B
	 * @param zB z coordinate of point B
	 * @return the Euclidean distance between points A and B
	 */
	double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
		return Math.sqrt(Math.pow(xB - xA, 2) + Math.pow(yB - yA, 2) + Math.pow(zB - zA, 2));
	}
	
	/**
	 * Total (kinetic + potential) energy of the system
	 * @return the total energy
	 */
	public double hamiltonian () {
		double energy = 0.0;
		for (int i = 0; i < np; i++) {
			Particle a = particles.get(i);
			energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass;
			for (int j = 0; j < np; j++) {
				if (i > j) {
					Particle b = particles.get(j);
					energy -= g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
				}
			}
		}
		return energy;
	}
	
	/**
	 * "Position outside, momentum between" variant
	 */
	public void solveQP () {
		integrator.solve(this, PhaseSpace.Q, PhaseSpace.P);
	}
	
	/**
	 * "Momentum outside, position between" variant (probably slightly slower, less accurate)
	 */
	public void solvePQ () {
		integrator.solve(this, PhaseSpace.P, PhaseSpace.Q);
	}
	
	/**
	 * Writes out current system state as JSON-formatted text
	 * @return the JSON string
	 */
	public String particlesJson () {
		StringBuilder json = new StringBuilder("[");
		Iterator<Particle> p = particles.iterator();
		json.append(p.next().toString());
		while (p.hasNext()) {
			json.append(',' + p.next().toString());
		}
		return (json + "]");
	}
}
