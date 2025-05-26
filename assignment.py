import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import art3d # Needed for 3D plotting
from scipy.integrate import simpson

class MobiusStrip:
    def __init__(self, R=1.0, w=0.3, n=200):
        
        # Input validation
        if not (isinstance(R, (int, float)) and R > 0):
            raise ValueError("Radius R must be a positive number.")
        if not (isinstance(w, (int, float)) and w > 0):
            raise ValueError("Width w must be a positive number.")
        if not (isinstance(n, int) and n > 1): # n > 1 to allow for gradients and diffs
            raise ValueError("Resolution n must be an integer greater than 1.")

        self.R = R
        self.w = w
        self.n = n

        # Generate the parameter space
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)

        # Compute the mesh grid (x, y, z)
        self.X, self.Y, self.Z = self.compute_surface()

    def compute_surface(self):
        
        U, V, R = self.U, self.V, self.R

        X = (R + V * np.cos(U / 2)) * np.cos(U)
        Y = (R + V * np.cos(U / 2)) * np.sin(U)
        Z = V * np.sin(U / 2)

        return X, Y, Z

    def surface_area(self):
        
        # Calculate step sizes for numerical differentiation
        du = self.u[1] - self.u[0]
        dv = self.v[1] - self.v[0]

        
        Xu, Xv = np.gradient(self.X, du, dv, edge_order=2)
        Yu, Yv = np.gradient(self.Y, du, dv, edge_order=2)
        Zu, Zv = np.gradient(self.Z, du, dv, edge_order=2)

        # Compute the components of the normal vector (cross product r_u x r_v)
        cross_x = Yu * Zv - Zu * Yv
        cross_y = Zu * Xv - Xu * Zv
        cross_z = Xu * Yv - Yu * Xv

        # Calculate the magnitude of the normal vector (dA element)
        dA_magnitude = np.sqrt(cross_x**2 + cross_y**2 + cross_z**2)

        
        area = simpson([simpson(row, self.v) for row in dA_magnitude], self.u)
        return area

    def edge_length(self):
        def length_of_curve(points):
            """Helper function to calculate the total length of a 3D curve."""
            # Calculate differences between consecutive points
            diffs = np.diff(points, axis=0)
            # Calculate Euclidean distance for each segment
            segment_lengths = np.linalg.norm(diffs, axis=1)
            # Sum all segment lengths
            return np.sum(segment_lengths)

        
        edge_points = np.array([self.X[-1, :], self.Y[-1, :], self.Z[-1, :]]).T

        # Calculate the length of this curve
        return length_of_curve(edge_points)

    def plot(self):
        """
        Plot the Möbius strip in 3D, including its surface and the single boundary edge.
        """
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the surface of the Mobius strip
        ax.plot_surface(self.X, self.Y, self.Z, rstride=1, cstride=1,
                        color='skyblue', edgecolor='none', alpha=0.8)

        
        edge_x = self.X[-1, :]
        edge_y = self.Y[-1, :]
        edge_z = self.Z[-1, :]
        ax.plot(edge_x, edge_y, edge_z, color='red', linewidth=3, label='Mobius Strip Edge')

        ax.set_title(f"Möbius Strip (R={self.R}, w={self.w}, n={self.n})")
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend()
        plt.tight_layout()
        plt.show()

# === Main Execution ===
if __name__ == "__main__":
    # Example usage with default parameters
    mobius = MobiusStrip(R=1.0, w=0.3, n=200)

    print(f"Approximate Surface Area: {mobius.surface_area():.4f}")
    print(f"Approximate Edge Length: {mobius.edge_length():.4f}")

    # Visualize the Möbius strip
    mobius.plot()

    # Example with different parameters
    print("\n--- Mobius Strip with R=2.0, w=0.5, n=100 ---")
    mobius2 = MobiusStrip(R=2.0, w=0.5, n=100)
    print(f"Approximate Surface Area: {mobius2.surface_area():.4f}")
    print(f"Approximate Edge Length: {mobius2.edge_length():.4f}")
    mobius2.plot()

    # Test edge cases or invalid inputs
    try:
        MobiusStrip(R=0, w=0.3, n=200)
    except ValueError as e:
        print(f"\nError creating MobiusStrip: {e}")

    try:
        MobiusStrip(R=1.0, w=-0.1, n=200)
    except ValueError as e:
        print(f"Error creating MobiusStrip: {e}")

    try:
        MobiusStrip(R=1.0, w=0.3, n=1)
    except ValueError as e:
        print(f"Error creating MobiusStrip: {e}")