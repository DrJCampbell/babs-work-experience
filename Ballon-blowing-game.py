import pygame, random, sys
from pygame.locals import *

# Initialisation
pygame.init()
pygame.font.init()

# Utility Functions
def dis_formula(a, b, c, d):
    dx = a - c
    dy = b - d
    return dx**2 + dy**2

# Display Setup
screen_w = 640
screen_h = 480
screen = pygame.display.set_mode((screen_w,screen_h))
pygame.display.set_caption("Balloon Blow-Up Game")
clock = pygame.time.Clock()

# Game Elements
balloon = pygame.image.load("Balloon.png").convert()
popped_balloon = pygame.image.load("Balloon_popped.png").convert()

width_b, height_b = balloon.get_size()
width_p_b, height_p_b = popped_balloon.get_size()
max_balloon_w = width_b * 1.5
max_balloon_h = height_b * 1.5

final_size = random.choice([max_balloon_w + 10, max_balloon_h + 10])
balloon_rect = balloon.get_rect(center=screen.get_rect().center)
og_width_b, og_height_b = balloon.get_size()
scale = min(screen_w / og_width_b , screen_h / og_height_b)
initial_scale = scale * 0.3
new_size = (int(og_width_b * initial_scale), int(og_height_b * initial_scale))

balloon = pygame.transform.smoothscale(balloon, new_size)
balloon_rect = balloon.get_rect(center=(screen_w // 2, screen_h // 2))
# Font & Text Setup
font = pygame.font.SysFont(None, 36)
balloon_scale = 1.01
# Game State
running = True
put_popped_balloon = False
# Main Game Loop
while running:
    # Event Handling
    for event in pygame.event.get():
        if event.type == QUIT:
            running = False
    if pygame.mouse.get_pressed()[0]:
        mouse_pos = pygame.mouse.get_pos()
        if balloon_rect.collidepoint(mouse_pos):
            balloon_scale *= 1.01  # Increase scale
            new_size = (int(og_width_b * balloon_scale), int(og_height_b * balloon_scale))
            balloon = pygame.transform.smoothscale(original_balloon, new_size)
            balloon_rect = balloon.get_rect(center=(screen_w // 2, screen_h // 2))  # Keep centered

            if new_size[0] > final_size or new_size[1] > final_size:
                put_popped_balloon = True

    # Drawing
    screen.fill((0,0,0))
    screen.blit(balloon, balloon_rect)
    if put_popped_balloon:
        screen.blit(popped_balloon, balloon_rect)

    #Update Display
    pygame.display.flip()
    clock.tick(60)
    
    
            
            
            
            
