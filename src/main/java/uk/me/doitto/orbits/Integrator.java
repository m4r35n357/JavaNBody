package uk.me.doitto.orbits;

/**
 * @author ian
 * 
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	EULER {
		@Override
		public void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 1.0);
			second.update(s, 1.0);
		}
	},
	STORMER_VERLET_2 {
		@Override
		public void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 0.5);
			second.update(s, 1.0);
			first.update(s, 0.5);
		}
	},
	STORMER_VERLET_4A {
		@Override
		public void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 0.675603595979828900000000);
			second.update(s, 1.351207191959657800000000);
			first.update(s, -0.175603595979828830000000);
			second.update(s, -1.702414383919315300000000);
			first.update(s, -0.175603595979828830000000);
			second.update(s, 1.351207191959657800000000);
			first.update(s, 0.675603595979828900000000);
		}
	},
	STORMER_VERLET_4 {
		@Override
		public void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			double c1, c2, c3, c4, d1, d2, d3;
			double cubeRoot2 = Math.pow(2.0, 1.0 / 3.0);
			double denom = 2.0 - cubeRoot2;
			c1 = c4 = 1.0 / (2.0 * denom);
			c2 = c3 = (1.0 - cubeRoot2) / (2.0 * denom);
			d1 = d3 = 1.0 / denom;
			d2 = - cubeRoot2 / denom;
			first.update(s, c1);
			second.update(s, d1);
			first.update(s, c2);
			second.update(s, d2);
			first.update(s, c3);
			second.update(s, d3);
			first.update(s, c4);
//			System.out.printf("c1: %.24f, c2: %.24f, d1: %.24f, d2: %.24f%n", c1, c2, d1, d2);
		}
	};

	public abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
