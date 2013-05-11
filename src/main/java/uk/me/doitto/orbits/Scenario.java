package uk.me.doitto.orbits;

import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_4;
import static uk.me.doitto.orbits.PhaseSpace.P;
import static uk.me.doitto.orbits.PhaseSpace.Q;

import java.util.ArrayList;
import java.util.List;

/**
 * @author ian
 * 
 * Scenarios
 */
public enum Scenario {
	CRISS_CROSS {
		@Override
		public void populate () {
			g = 0.05;
			ts = 0.001;
			outputInterval = 1000;
			simulationTime = 1.0e4;
			bodies.add(new Particle(1.07590, 0.0, 0.0, 0.1, 0.1, 0.0, 1.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		}
	},
	TWO_BODY {
		@Override
		public void populate () {
			g = 0.05;
			ts = 0.001;
			outputInterval = 1000;
			simulationTime = 1.0e5;
			bodies.add(new Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		}
	},
	THREE_BODY {
		@Override
		public void populate () {
			g = 1.0;
			ts = 0.001;
			outputInterval = 1000;
			simulationTime = 1.0e5;
			bodies.add(new Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0));
			bodies.add(new Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0));
			bodies.add(new Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0));
		}
	},
	FOUR_BODY {
		@Override
		public void populate () {
			g = 3.5;
			ts = 0.001;
			outputInterval = 1000;
			simulationTime = 1.0e5;
			bodies.add(new Particle(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0));
			bodies.add(new Particle(-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0));
			bodies.add(new Particle(1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0));
			bodies.add(new Particle(-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0));
		}
	},
	EIGHT_BODY {
		@Override
		public void populate () {
			g = 0.05;
			ts = 0.001;
			outputInterval = 1000;
			simulationTime = 1.0e4;
			bodies.add(new Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0));
			bodies.add(new Particle(0.0, 4.5, 0.4, -0.2, 0.0, 1.8, 2.0));
			bodies.add(new Particle(-6.0, 0.0, -0.4, 0.0, -0.6, 1.0, 3.0));
			bodies.add(new Particle(3.0, 0.0, -0.2, 0.0, 5.8, -0.2, 5.0));
			bodies.add(new Particle(0.0, -4.0, 0.1, -3.6, 0.0, 0.2, 4.0));
			bodies.add(new Particle(-4.0, 0.0, -0.1, 0.0, -0.2, -2.6, 3.0));
			bodies.add(new Particle(8.0, 0.0, -0.3, 0.0, 1.2, -0.2, 3.0));
			bodies.add(new Particle(0.0, 4.0, -0.2, -4.8, 0.0, -0.2, 4.0));
		}
	};
	
	double g = 0.0;
	
	double simulationTime = 0.0;

	double ts = 0.0;
	
	int outputInterval;
	
	List<Particle> bodies = new ArrayList<Particle>();

	private Scenario () {
		populate();
	}
	
	public abstract void populate ();

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
			if ((n % s.outputInterval) == 0) {
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
