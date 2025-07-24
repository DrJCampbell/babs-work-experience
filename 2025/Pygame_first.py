import pygame

pygame.init()
# Set up the display
width, height = 800, 600
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Balloon Game")
clock = pygame.time.Clock()
# Load images
balloon_image = pygame.image.load("C:\\Users\\workexperience25\\Documents\\GitHub\\babs-work-experience\\Balloon.png")
balloon_popped_image = pygame.image.load("C:\\Users\\workexperience25\\Documents\\GitHub\\babs-work-experience\\Balloon_popped.png")

# Set original balloon size
balloon_width, balloon_height = 200, 200  # Set your desired starting size here

# Centre the balloon
balloon_x = (width - balloon_width) // 2
balloon_y = (height - balloon_height) // 2

max_width = balloon_width * 2  # Max size before popping
max_height = balloon_height * 2
popped = False

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))

    mouse_pos = pygame.mouse.get_pos()
    mouse_pressed = pygame.mouse.get_pressed()
    balloon_rect = pygame.Rect(balloon_x, balloon_y, balloon_width, balloon_height)

    if not popped:
        # Check if mouse is pressed on the balloon
        if balloon_rect.collidepoint(mouse_pos) and mouse_pressed[0]:
            balloon_width += 1
            balloon_height += 1
            # Re-center as the balloon grows
            balloon_x = (width - balloon_width) // 2
            balloon_y = (height - balloon_height) // 2
            if balloon_width > max_width or balloon_height > max_height:
                popped = True

        # Draw the (possibly scaled) balloon
        scaled_balloon = pygame.transform.scale(balloon_image, (balloon_width, balloon_height))
        screen.blit(scaled_balloon, (balloon_x, balloon_y))
    else:
        # Draw the popped balloon image, centered
        popped_x = (width - balloon_popped_image.get_width()) // 2
        popped_y = (height - balloon_popped_image.get_height()) // 2
        screen.blit(balloon_popped_image, (popped_x, popped_y))