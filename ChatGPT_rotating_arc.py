import tkinter as tk
import math

class RotatingArcApp:
    def __init__(self, root):
        self.root = root
        self.canvas = tk.Canvas(root, width=400, height=400, bg="white")
        self.canvas.pack()

        # Center and initial dimensions of the arc
        self.center_x = 200
        self.center_y = 200
        self.radius = 50

        # Initial angle for rotation
        self.angle = 0

        # Create the initial arc
        self.arc = self.create_arc(self.angle)

        # Start rotating
        self.rotate_arc()

    def create_arc(self, angle):
        """Create an arc at a specific angle."""
        # Convert angle to radians
        rad_angle = math.radians(angle)

        # Calculate the new bounding box based on rotation
        x0 = self.center_x - self.radius
        y0 = self.center_y - self.radius
        x1 = self.center_x + self.radius
        y1 = self.center_y + self.radius

        # Adjust start and extent angles for rotation
        start_angle = (angle + 0) % 360  # Adjust as needed
        extent_angle = 90  # Arc extent (e.g., quarter-circle)

        # Create an arc
        return self.canvas.create_arc(x0, y0, x1, y1, start=start_angle, extent=extent_angle, style=tk.ARC, outline="blue", width=2)

    def rotate_arc(self):
        """Rotate the arc incrementally."""
        # Delete the old arc
        self.canvas.delete(self.arc)

        # Increment the angle
        self.angle = (self.angle + 5) % 360  # Rotate by 5 degrees each frame

        # Create a new arc with the updated angle
        self.arc = self.create_arc(self.angle)

        # Schedule the next rotation
        self.root.after(50, self.rotate_arc)  # Rotate every 50ms

# Set up the main application window
root = tk.Tk()
root.title("Rotating Arc")
app = RotatingArcApp(root)
root.mainloop()
