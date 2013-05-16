package uk.me.doitto.orbits;

/**
 * @author ian
 * 
 * Definitions of the two Hamiltonian coordinates, position and momentum, specifically their update methods
 */
public enum PhaseSpace {	
	Q {  // Position
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
	P1 {  // Momentum
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				double am = a.mass;
				double aQx = a.qX;
				double aQy = a.qY;
				double aQz = a.qZ;
				for (int j = 0; j < s.np; j++) {
					Particle b = s.particles.get(j);
					double bQx = b.qX;
					double bQy = b.qY;
					double bQz = b.qZ;
					if (i > j) {
						double tmp = - c * s.g * am * b.mass / Math.pow(s.distance(aQx, aQy, aQz, bQx, bQy, bQz), 3) * s.timeStep;
						double dPx = (bQx - aQx) * tmp;
						double dPy = (bQy - aQy) * tmp;
						double dPz = (bQz - aQz) * tmp;
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
	},
	P {  // Momentum
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				for (int j = 0; j < s.np; j++) {
					Particle b = s.particles.get(j);
					if (i > j) {
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

	abstract void update (Symplectic symplectic, double coefficient);
}
