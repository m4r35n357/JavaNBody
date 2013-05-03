package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
public enum Coordinates implements Coordinate {
	POSITION {
		@Override
		public void update (Symplectic symplectic, double coefficient) {
			Particle a;
			int i;
			double tmp;
			for (i = 0; i < symplectic.np; i++) {
				a = symplectic.particles.get(i);
				tmp = coefficient / a.mass * symplectic.ts;
				a.qX += a.pX * tmp;
				a.qY += a.pY * tmp;
				a.qZ += a.pZ * tmp;
			}
		}
	},
	MOMENTUM {
		@Override
		public void update (Symplectic symplectic, double coefficient) {
			Particle a, b;
			int i, j;
			double tmp, dPx, dPy, dPz;
			for (i = 0; i < symplectic.np; i++) {
				a = symplectic.particles.get(i);
				for (j = 0; j < symplectic.np; j++) {
					b = symplectic.particles.get(j);
					tmp = - coefficient * symplectic.g * a.mass * b.mass / Math.pow(symplectic.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * symplectic.ts;
					if (i > j) {
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

	@Override
	public abstract void update (Symplectic symplectic, double coefficient);
}
