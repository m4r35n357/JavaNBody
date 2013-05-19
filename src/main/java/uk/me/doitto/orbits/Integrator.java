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
		}
	},
	STORMER_VERLET_2 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2 (s, first, second, 1.0);
		}
	},
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 1.351207191959657);
			symplectic2(s, first, second, -1.702414383919315);
			symplectic2(s, first, second, 1.351207191959657);
		}
	},
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.784513610477560e0);
			symplectic2(s, first, second, 0.235573213359357e0);
			symplectic2(s, first, second, -1.17767998417887e0);
			symplectic2(s, first, second, 1.31518632068391e0);
			symplectic2(s, first, second, -1.17767998417887e0);
			symplectic2(s, first, second, 0.235573213359357e0);
			symplectic2(s, first, second, 0.784513610477560e0);
		}
	},
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.104242620869991e1);
			symplectic2(s, first, second, 0.182020630970714e1);
			symplectic2(s, first, second, 0.157739928123617e0);
			symplectic2(s, first, second, 0.244002732616735e1);
			symplectic2(s, first, second, -0.716989419708120e-2);
			symplectic2(s, first, second, -0.244699182370524e1);
			symplectic2(s, first, second, -0.161582374150097e1);
			symplectic2(s, first, second, -0.17808286265894516e1);
			symplectic2(s, first, second, -0.161582374150097e1);
			symplectic2(s, first, second, -0.244699182370524e1);
			symplectic2(s, first, second, -0.716989419708120e-2);
			symplectic2(s, first, second, 0.244002732616735e1);
			symplectic2(s, first, second, 0.157739928123617e0);
			symplectic2(s, first, second, 0.182020630970714e1);
			symplectic2(s, first, second, 0.104242620869991e1);
		}
	};

	protected final void symplectic2 (Symplectic s, PhaseSpace first, PhaseSpace second, double step) {
		first.update(s, 0.5 * step);
		second.update(s, step);
		first.update(s, 0.5 * step);
	}
	
	abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
