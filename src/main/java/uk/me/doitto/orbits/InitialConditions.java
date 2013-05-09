package uk.me.doitto.orbits;

import java.util.ArrayList;
import java.util.List;

/**
 * @author ian
 * 
 * Scenarios
 */
public enum InitialConditions {
	CRISS_CROSS {
		@Override
		public void populate () {
			g = 0.05;
			ts = 0.001;
			simulationTime = 10000.0;
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
			simulationTime = 10000.0;
			bodies.add(new Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		}
	},
	THREE_BODY {
		@Override
		public void populate () {
			g = 1.0;
			ts = 0.001;
			simulationTime = 1.0e6;
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
			simulationTime = 10000.0;
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
	
	List<Particle> bodies = new ArrayList<Particle>();

	private InitialConditions () {
		populate();
	}
	
	public abstract void populate ();
}
