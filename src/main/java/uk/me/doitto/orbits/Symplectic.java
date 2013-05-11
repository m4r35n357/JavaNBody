package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.List;

/**
 * @author ian
 *
 * Top level class for symplectic integrator simulations
 */
public class Symplectic {
	
	double iterations;
	
	double g;
	
	double timeStep;
	
	int outputInterval;
	
	int np = 0;
	
	List<Particle> particles;
	
	Integrator integrator;
	
	public Symplectic (Scenario ic, Integrator integrator) {
		this.particles = ic.bodies;
		this.np = ic.bodies.size();
		this.g = ic.g;
		this.timeStep = ic.ts;
		this.outputInterval = ic.outputInterval;
		this.iterations = ic.simulationTime / ic.ts;
		this.integrator = integrator;
	}
	
	public Symplectic (double g, double simulationTime, double timeStep, int outputInterval, List<Particle> bodies, String integrator) {
		// destroy reference to input array so the client can't change i
		this.particles = new ArrayList<Particle>(bodies);
		this.np = bodies.size();
		this.g = g;
		this.timeStep = timeStep;
		this.outputInterval = outputInterval;
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
}
