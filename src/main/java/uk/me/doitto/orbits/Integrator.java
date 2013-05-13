package uk.me.doitto.orbits;

/**
 * @author ian
 * 
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	EULER {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 1.0);
			second.update(s, 1.0);
			s.cog();
		}
	},
	STORMER_VERLET_2 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 0.5);
			second.update(s, 1.0);
			first.update(s, 0.5);
			s.cog();
		}
	},
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 0.675603595979828900000000);
			second.update(s, 1.351207191959657800000000);
			first.update(s, -0.175603595979828830000000);
			second.update(s, -1.702414383919315300000000);
			first.update(s, -0.175603595979828830000000);
			second.update(s, 1.351207191959657800000000);
			first.update(s, 0.675603595979828900000000);
			s.cog();
		}
	};

	abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
