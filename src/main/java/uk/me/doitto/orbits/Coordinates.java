package uk.me.doitto.orbits;

/**
 * @author ian
 * 
 * Definitions of the two Hamiltonian coordinates, position and momentum, specifically their update methods
 */
public enum Coordinates {
	// Position
	Q {
		@Override
		public void update (Symplectic s, double c) {
			Particle a;
			double tmp;
			for (int i = 0; i < s.np; i++) {
				a = s.particles.get(i);
				tmp = c / a.mass * s.timeStep;
				a.qX += a.pX * tmp;
				a.qY += a.pY * tmp;
				a.qZ += a.pZ * tmp;
			}
		}
	},
	// Momentum
	P {
		@Override
		public void update (Symplectic s, double c) {
			Particle a, b;
			double tmp, dPx, dPy, dPz;
			for (int i = 0; i < s.np; i++) {
				a = s.particles.get(i);
				for (int j = 0; j < s.np; j++) {
					b = s.particles.get(j);
					if (i > j) {
						tmp = - c * s.g * a.mass * b.mass / Math.pow(s.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * s.timeStep;
						dPx = (b.qX - a.qX) * tmp;
						dPy = (b.qY - a.qY) * tmp;
						dPz = (b.qZ - a.qZ) * tmp;
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

	public abstract void update (Symplectic symplectic, double coefficient);
}
