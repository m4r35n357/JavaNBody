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
public enum InitialConditions implements IInitialConditions {
	CRISS_CROSS {
		@Override
		public void populate () {
			bodies.add(new Particle(1.07590, 0.0, 0.0, 0.1, 0.1, 0.0, 1.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		}
	},
	TWO_BODY {
		@Override
		public void populate () {
			bodies.add(new Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0));
			bodies.add(new Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0));
		}
	},
	EIGHT_BODY {
		@Override
		public void populate () {
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
	
	public List<Particle> bodies = new ArrayList<Particle>();

	private InitialConditions () {
		populate();
	}
	
	@Override
	public abstract void populate ();
}
