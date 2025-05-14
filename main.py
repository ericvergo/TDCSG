from math import pi, cos, sin, hypot, sqrt, floor
from PIL import Image, ImageDraw
import colorsys
from time import sleep


class PortraitGenerator:
    def __init__(self, n, alpha, beta, radius2, disc_distance, resolution):
        """Initialize the portrait generator with geometric parameters.

        Args:
            n: Base value for rotation calculations
            alpha: Control parameter for left rotation
            beta: Control parameter for right rotation
            radius2: Radius of the second disc
            disc_distance: Distance between disc centers
            resolution: Pixels per unit distance (determines image quality)
        """
        # Core geometric parameters
        self.n = n
        self.alpha = alpha
        self.beta = beta
        self.radius2 = radius2
        self.disc_distance = disc_distance
        self.resolution = resolution

        # Derived values
        self.pixel_shift = 0.5 / resolution
        self.multiplier = 10 ** 11  # For float comparison precision

        # Complex number rotations
        angle = 2 * pi / n
        self.left_rotation = complex(cos(angle), sin(angle)) ** alpha
        self.right_rotation = complex(cos(-angle), sin(-angle)) ** beta
        self.shift_complex = complex(disc_distance, 0)

        # Image parameters
        self.x_shift = resolution + 2
        self.y_shift = resolution + 2
        x_size = floor((1 + radius2 + disc_distance) * resolution) + 4
        y_size = floor(2 * resolution) + 4

        # Initialize image and tracking structures
        self.image = Image.new(mode="RGB", size=(int(x_size), int(y_size)), color=(255, 255, 255))
        self.draw = ImageDraw.Draw(self.image)
        self.visited_points = {}
        self.points_to_process = []
        self.fill_points = []
        self.iteration_count = 0

    def left_transform(self, point):
        """Rotate point around the origin (first disc center)."""
        if self.rotation_type == "simple generator":
            return point * self.left_rotation
        else:
            result = []
            for i in range(floor(self.n / self.alpha)):
                point = point * self.left_rotation
                result.append(point)
            return result

    def right_transform(self, point):
        """Rotate point around the second disc center."""
        shifted = point - self.shift_complex
        if self.rotation_type == "simple generator":
            return shifted * self.right_rotation + self.shift_complex
        else:
            result = []
            for i in range(floor(self.n / self.beta)):
                shifted = shifted * self.right_rotation
                result.append(shifted + self.shift_complex)
            return result

    def add_point(self, points, mode):
        """Add points to processing queue and optionally draw them.

        Args:
            points: List of points to process
            mode: 0 for drawing points, 1 for collecting fill points
        """
        for point in points:
            # Create a unique identifier for this point
            point_id = f"{round(point.real * self.multiplier)}{round(point.imag * self.multiplier)}"

            # Skip if we've already processed this point
            if point_id in self.visited_points:
                continue

            self.visited_points[point_id] = 0

            # Handle different rotation types
            if self.rotation_type == "simple generator":
                self.points_to_process.append(point)
            else:
                # Check if point is in the intersection of the discs
                if (hypot(point.real, point.imag) <= 1 and
                        hypot(point.real - self.disc_distance, point.imag) <= self.radius2):
                    self.points_to_process.append(point)

            # Convert complex coordinates to pixel coordinates
            x = int(point.real * self.resolution) + self.x_shift
            y = int(point.imag * self.resolution) + self.y_shift

            if mode == 0:  # Draw point on the boundary
                self.draw.point((x, y), fill=(0, 0, 0))
            elif mode == 1:  # Add to fill list for coloring
                self.fill_points.append((x, y))

    def process_points(self, mode):
        """Process all points in the queue, applying transformations.

        Args:
            mode: 0 for drawing points, 1 for collecting fill points
        """
        while self.points_to_process and self.iteration_count < self.max_iterations:
            current_point = self.points_to_process.pop(0)

            if self.rotation_type == "simple generator":
                # Apply appropriate transformation based on point location
                if hypot(current_point.real, current_point.imag) <= 1:
                    current_point = self.left_transform(current_point)
                if hypot(current_point.real - self.disc_distance, current_point.imag) <= self.radius2:
                    current_point = self.right_transform(current_point)
                self.add_point([current_point], mode)
            else:
                # Apply both transformations and add resulting points
                self.add_point(self.left_transform(current_point), mode)
                self.add_point(self.right_transform(current_point), mode)

            self.iteration_count += 1

    def fill_regions(self):
        """Find and color regions based on iteration counts."""
        # Define a neighborhood of pixels to check
        neighbors = [(x, y) for x in range(-3, 4) for y in range(-3, 4)]

        # Check pixels in the intersection region
        for x in range(int(self.resolution * (1 + self.disc_distance - self.radius2)), 2 * self.resolution, 1):
            for y in range(0, self.image.height, 1):
                # Check if pixel is in the intersection of both discs
                x_coord = x / self.resolution - 1 + self.pixel_shift
                y_coord = y / self.resolution - 1 - self.pixel_shift

                if (hypot(x_coord, y_coord) <= 1 and
                        hypot(x_coord - self.disc_distance, y_coord) <= self.radius2):

                    # Skip if pixel is already colored
                    if self.image.getpixel((x, y)) != (255, 255, 255):
                        continue

                    # Check if this is an isolated white pixel
                    is_isolated = True
                    for dx, dy in neighbors:
                        try:
                            if self.image.getpixel((x + dx, y + dy)) != (255, 255, 255):
                                is_isolated = False
                                break
                        except IndexError:
                            continue

                    if is_isolated:
                        self.fill_points = [(x, y)]
                        self.points_to_process.append(complex(x_coord, y_coord))
                        self.iteration_count = 0
                        self.process_points(1)

                        # Color based on iteration count
                        hue = 0.66 * (sqrt(self.iteration_count * 700) -
                                      floor(sqrt(self.iteration_count * 700)))
                        r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 0.75)
                        fill_color = (int(255 * r), int(225 * g), int(255 * b))

                        for fill_point in self.fill_points:
                            ImageDraw.floodfill(self.image, fill_point, fill_color, (0,0,0))

    def generate(self, point_set, rotation_type, fill=True):
        """Generate the portrait image.

        Args:
            point_set: "full boundary" or a specific starting point as [x, y]
            rotation_type: "simple generator" or "full unbandaging"
            fill: Whether to color the regions
        """
        self.rotation_type = rotation_type
        self.points_to_process = []
        self.fill_points = []
        self.visited_points = {}
        self.iteration_count = 0

        # Set maximum iterations based on point set
        if point_set == "full boundary":
            self.max_iterations = 10 ** 20
            # Sample points along the boundaries of both discs
            max_samples = max(2, 2 * self.n) * self.resolution
            for i in range(int(2 * max_samples)):
                angle = pi * i / max_samples
                # Points on the first disc boundary
                x, y = cos(angle) * 0.9999999999, sin(angle) * 0.9999999999
                self.add_point([complex(x, y)], 0)
                # Points on the second disc boundary
                x, y = x * self.radius2, y * self.radius2
                self.add_point([complex(x + self.disc_distance, y)], 0)
        else:
            # Start with a single specified point
            self.add_point([complex(point_set[0], point_set[1])], 0)
            self.max_iterations = (self.resolution * 2.5) ** 2

        # Process boundary points
        self.process_points(0)

        # Fill regions if requested
        if fill:
            self.fill_regions()

        # Display the image
        self.image.show()

        return self.image


# Main execution
if __name__ == "__main__":
    # Distance between disc centers
    distance = 1.25

    # Resolution (pixels per unit distance)
    resolution = 1000

    # Create first portrait with simple generator rotation
    generator1 = PortraitGenerator(
        n=8,
        alpha=1,
        beta=3,
        radius2=1,
        disc_distance=distance,
        resolution=resolution
    )
    generator1.generate(point_set="full boundary", rotation_type="simple generator", fill=True)


    # Distance between disc centers
    distance = 0.5

    # Resolution (pixels per unit distance)
    resolution = 500

    # Create second portrait with full unbandaging rotation
    generator2 = PortraitGenerator(
        n=6,
        alpha=1,
        beta=1,
        radius2=1,
        disc_distance=distance,
        resolution=resolution
    )
    generator2.generate(point_set="full boundary", rotation_type="full unbandaging", fill=True)

    # Distance between disc centers
    distance = sqrt(2)

    # Resolution (pixels per unit distance)
    resolution = 1000

    # Create third portrait with a single point under a simple generator
    generator2 = PortraitGenerator(
        n=12,
        alpha=1,
        beta=1,
        radius2=1,
        disc_distance=distance,
        resolution=resolution
    )
    generator2.generate(point_set=(distance/2, 0), rotation_type="simple generator", fill=False)