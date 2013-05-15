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
			this.stormerVerletBase (s, first, second, 1.0);
//			first.update(s, 0.5);
//			second.update(s, 1.0);
//			first.update(s, 0.5);
			s.cog();
		}
	},
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			this.stormerVerletBase(s, first, second, 1.351207191959657);
			this.stormerVerletBase(s, first, second, -1.702414383919315);
			this.stormerVerletBase(s, first, second, 1.351207191959657);
//			first.update(s, 0.675603595979828900000000);
//			second.update(s, 1.351207191959657800000000);
//			first.update(s, -0.175603595979828830000000);
//			second.update(s, -1.702414383919315300000000);
//			first.update(s, -0.175603595979828830000000);
//			second.update(s, 1.351207191959657800000000);
//			first.update(s, 0.675603595979828900000000);
			s.cog();
		}
	},
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			this.stormerVerletBase(s, first, second, 0.784513610477560e0);
			this.stormerVerletBase(s, first, second, 0.235573213359357e0);
			this.stormerVerletBase(s, first, second, -1.17767998417887e0);
			this.stormerVerletBase(s, first, second, 1.31518632068391e0);
			this.stormerVerletBase(s, first, second, -1.17767998417887e0);
			this.stormerVerletBase(s, first, second, 0.235573213359357e0);
			this.stormerVerletBase(s, first, second, 0.784513610477560e0);
			s.cog();
		}
	},
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			this.stormerVerletBase(s, first, second, 0.104242620869991e1);
			this.stormerVerletBase(s, first, second, 0.182020630970714e1);
			this.stormerVerletBase(s, first, second, 0.157739928123617e0);
			this.stormerVerletBase(s, first, second, 0.244002732616735e1);
			this.stormerVerletBase(s, first, second, -0.716989419708120e-2);
			this.stormerVerletBase(s, first, second, -0.244699182370524e1);
			this.stormerVerletBase(s, first, second, -0.161582374150097e1);
			this.stormerVerletBase(s, first, second, -0.17808286265894516e1);
			this.stormerVerletBase(s, first, second, -0.161582374150097e1);
			this.stormerVerletBase(s, first, second, -0.244699182370524e1);
			this.stormerVerletBase(s, first, second, -0.716989419708120e-2);
			this.stormerVerletBase(s, first, second, 0.244002732616735e1);
			this.stormerVerletBase(s, first, second, 0.157739928123617e0);
			this.stormerVerletBase(s, first, second, 0.182020630970714e1);
			this.stormerVerletBase(s, first, second, 0.104242620869991e1);
			s.cog();
		}
	};

	protected final void stormerVerletBase (Symplectic s, PhaseSpace first, PhaseSpace second, double step) {
		first.update(s, 0.5 * step);
		second.update(s, step);
		first.update(s, 0.5 * step);
	}
	
	abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
