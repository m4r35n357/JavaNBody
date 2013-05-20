package uk.me.doitto.orbits;

/**
 * @author ian
 * <p>
 * Definitions of the two Hamiltonian coordinates; position and momentum, specifically their update methods
 */
public enum PhaseSpace {	
	/**
	 * Position
	 * <p>
	 * update implements dH/dp, which in this case is a function of p only
	 */
	Q {
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				double tmp = c / a.mass * s.timeStep;
				a.qX += a.pX * tmp;
				a.qY += a.pY * tmp;
				a.qZ += a.pZ * tmp;
			}
		}
	},	
	/**
	 * Momentum, update implements -dH/dq, which in this case is a function of q only
	 */
	P {
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				for (int j = 0; j < s.np; j++) {
					if (i > j) {
						Particle b = s.particles.get(j);
						double tmp = - c * s.g * a.mass * b.mass / Math.pow(s.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * s.timeStep;
						double dPx = (b.qX - a.qX) * tmp;
						double dPy = (b.qY - a.qY) * tmp;
						double dPz = (b.qZ - a.qZ) * tmp;
						a.pX -= dPx;
						a.pY -= dPy;
						a.pZ -= dPz;
						b.pX += dPx;
						b.pY += dPy;
						b.pZ += dPz;
					}
				}
			}
		}
	};

	/**
	 * Update one component of the system state; either position or momentum
	 * @param s the Symplectic instance
	 * @param coefficient the Yoshida numbers for composed Stormer-Verlet methods
	 */
	abstract void update (Symplectic s, double coefficient);
}
