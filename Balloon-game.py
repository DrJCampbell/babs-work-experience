import pygame, sys, random
from pygame.locals import *

# Initialization
pygame.init()
pygame.font.init()

# Global colour library
colors = {
    "black":   (0, 0, 0),
    "white":   (255, 255, 255),
    "red":     (255, 0, 0),
    "green":   (0, 255, 0),
    "blue":    (0, 0, 255),
    "yellow":  (255, 255, 0),
    "cyan":    (0, 255, 255),
    "magenta": (255, 0, 255),
    "gray":    (128, 128, 128),
    "orange":  (255, 165, 0),
    "purple":  (128, 0, 128),
    "brown":   (165, 42, 42),
    "pink":    (255, 192, 203),
    "midnight_blue": (25, 25, 112)
}

# Utility Functions
def dis_formula(a, b, c, d):
    dx = a - c
    dy = b - d
    return dx**2 + dy**2

def ran_color():
    ran_choice = random.choice(list(colors.values()))
    while ran_choice == colors["midnight_blue"]:
        ran_choice = ran_color()
    return ran_choice

# Display Setup
screen_width = 640
screen_height = 480
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption("Balloon Colour Game")
clock = pygame.time.Clock()

# Game Elements
radius = 75
x = screen_width // 2
y = screen_height // 2
current_color = ran_color()

# Font & Text Setup
font = pygame.font.SysFont(None, 36)

# Game State
score = 0
highscore = 0

#  Main Game Loop
while True:
    # Event Handling
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit(0)
        if event.type == MOUSEBUTTONDOWN:
            mouse_x, mouse_y = event.pos
            if dis_formula(mouse_x, mouse_y, x, y) <= radius**2:
                current_color = ran_color()
                if current_color != colors["black"]:
                    score += 1
                else:
                    if score > highscore:
                        highscore = score
                    score = 0

    # Drawing
    screen.fill(colors["midnight_blue"])
    pygame.draw.circle(screen, current_color, (x, y), radius)

    # Score Display
    score_text = font.render(f"Score: {score}", True, colors["white"])
    score_rect = score_text.get_rect(topright=(screen_width - 10, 10))
    screen.blit(score_text, score_rect)

    highscore_text = font.render(f"Highscore: {highscore}", True, colors["white"])
    highscore_rect = highscore_text.get_rect(topleft=(10, 50))
    screen.blit(highscore_text, highscore_rect)

    # Update Display
    pygame.display.update()
    clock.tick(60)





