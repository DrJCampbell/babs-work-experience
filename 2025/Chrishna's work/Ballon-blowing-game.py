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
popped_balloon = pygame.image.load("Balloon_popped.png").convert_alpha()
original_balloon = balloon.copy()

width_b, height_b = balloon.get_size()
print(width_b, height_b)
width_p_b, height_p_b = popped_balloon.get_size()
print(width_p_b, height_p_b)
max_balloon_w = width_b * 0.22
max_balloon_h = height_b * 0.22

final_size_w = random.choice([max_balloon_w + 10, max_balloon_h + 20])
final_size_h = random.choice([max_balloon_w + 10, max_balloon_h + 20])
final_size_dimen = (final_size_w, final_size_h)
balloon_rect = balloon.get_rect(center=screen.get_rect().center)
og_width_b, og_height_b = original_balloon.get_size()
scale = min(screen_w / og_width_b , screen_h / og_height_b)
initial_scale = scale * 0.3
balloon_scale = initial_scale

new_size = (int(og_width_b * balloon_scale), int(og_height_b * balloon_scale))
balloon = pygame.transform.smoothscale(original_balloon, new_size)
balloon_rect = balloon.get_rect(center=(screen_w // 2, screen_h // 2))
# Font & Text Setup
font = pygame.font.SysFont(None, 36)

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
            balloon_scale += 0.01  # Increase scale
            new_size = (int(og_width_b * balloon_scale), int(og_height_b * balloon_scale))
            balloon = pygame.transform.smoothscale(original_balloon, new_size)
            balloon_rect = balloon.get_rect(center=(screen_w // 2, screen_h // 2))  # Keep centered
            print(balloon.get_size()[0])
            print(final_size_dimen[0])
            if balloon.get_size()[0] > max_balloon_w or balloon.get_size()[1] > max_balloon_h:
                put_popped_balloon = True
    # Drawing
    screen.fill((0,0,0))
    
    if put_popped_balloon:
    # Scale popped image to match balloon size
        popped_scaled = pygame.transform.smoothscale(popped_balloon, balloon.get_size())
        popped_rect = popped_scaled.get_rect(center=(screen_w // 2, screen_h // 2))
        screen.blit(popped_scaled, popped_rect)
    else:
        screen.blit(balloon, balloon_rect)

    #Update Display
    pygame.display.flip()
    clock.tick(60)
    
    
            
            
            
            
