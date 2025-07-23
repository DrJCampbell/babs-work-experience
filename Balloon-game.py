import pygame, sys, random
from pygame.locals import *
pygame.init()
pygame.font.init()

#Global dictionary
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
    "pink":    (255, 192, 203)
}


def dis_formula(a,b,c,d):
    dx = a - c
    dy = b - d
    distance_squared = dx**2 + dy**2
    return distance_squared

def ran_color(): # To be used with the colors dictionary
    ran_choice = random.choice(list(colors.values()))
    return ran_choice
    

screen_width = 640
screen_height = 480
radius = 75

font = pygame.font.SysFont(None, 36)
leaderboard = []
score = 0

screen = pygame.display.set_mode((screen_width,
                                  screen_height))
pygame.display.set_caption("Balloon luck game")
screen.fill((0,0,0))
current_color = ran_color()


x = screen_width // 2
y = screen_height // 2


while True:

    for event in pygame.event.get():
        if (event.type == QUIT):
            pygame.quit()
            sys.exit(0)
        if (event.type == MOUSEBUTTONDOWN):
            mouse_x, mouse_y = event.pos
            if dis_formula(mouse_x, mouse_y, x, y) <= radius**2:
                current_color = ran_color()
                while current_color == colors["midnight_blue"]:
                    current_color = ran_color()
                    

                if current_color != colors["black"]:
                    score += 1
                    
                else:
                    leaderboard.append(score)
                    score = 0

            else:
                pass
    colors["midnight_blue"] = (25, 25, 112)
    screen.fill(colors["midnight_blue"])

    pygame.draw.circle(screen, current_color, (x,y), radius)
    score_text = font.render(f"Score: {score}", True, colors["white"])
    text_rect = score_text.get_rect(topright = (screen_width - 10, 10))
    if len(leaderboard) == 0:
        i = " "
    else:
        i = leaderboard[-1]

    
    score_text2 = font.render(f"Highscore: {i}", True, colors["white"])
    text_rect2 = score_text2.get_rect(topleft = (10,50))

    screen.blit(score_text, text_rect)
    screen.blit(score_text2, text_rect2)
                
    pygame.display.update()
    

    

