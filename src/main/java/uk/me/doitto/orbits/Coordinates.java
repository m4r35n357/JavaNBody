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
			int i;
			double tmp;
			for (i = 0; i < s.np; i++) {
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
			double separation;
			int i, j;
			double tmp, dPx, dPy, dPz;
			for (i = 0; i < s.np; i++) {
				a = s.particles.get(i);
				for (j = 0; j < s.np; j++) {
					b = s.particles.get(j);
					if (i > j) {
						separation = s.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
						if (separation > (a.radius() + b.radius())) {
							tmp = - c * s.g * a.mass * b.mass / Math.pow(separation, 3) * s.timeStep;
						} else {
							tmp = - c * s.g * a.mass * b.mass / Math.pow(a.radius() + b.radius(), 3) * s.timeStep;
						}
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
