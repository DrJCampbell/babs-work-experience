# Example file showing a circle moving on screen
import pygame # type: ignore
import sys

# pygame setup
pygame.init()
screen = pygame.display.set_mode((1280, 720))
clock = pygame.time.Clock()
running = True
dt = 0

player1_pos = pygame.Vector2(screen.get_width() / 2, screen.get_height() / 2)
player2_pos = pygame.Vector2(screen.get_width() / 2, screen.get_height() / 2)




while running:
    # poll for events
    # pygame.QUIT event means the user clicked X to close your window
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # fill the screen with a color to wipe away anything from last frame
    screen.fill("blue")

    pygame.draw.circle(screen, "red", player1_pos, 10)

    keys = pygame.key.get_pressed()
    if keys[pygame.K_w]:
        player1_pos.y -= 300 * dt
    if keys[pygame.K_s]:
        player1_pos.y += 300 * dt
    if keys[pygame.K_a]:
        player1_pos.x -= 300 * dt
    if keys[pygame.K_d]:
        player1_pos.x += 300 * dt
    
    pygame.draw.circle(screen, "purple", player2_pos, 10)
   
    keys = pygame.key.get_pressed()
    if keys[pygame.K_UP]:
        player2_pos.y -= 300 * dt
    if keys[pygame.K_DOWN]:
        player2_pos.y += 300 * dt
    if keys[pygame.K_LEFT]:
        player2_pos.x -= 300 * dt
    if keys[pygame.K_RIGHT]:
        player2_pos.x += 300 * dt

    ball_radius = 10
    speed = 5

    player1_pos = pygame.Vector2(600, 300)
    player2_pos = pygame.Vector2(600, 300)
   
    player1_pos[0] = max(ball_radius, min(1280 - ball_radius, player1_pos[0]))
    player1_pos[1] = max(ball_radius, min(720 - ball_radius, player1_pos[1]))
    player2_pos[0] = max(ball_radius, min(1280 - ball_radius, player2_pos[0]))
    player2_pos[1] = max(ball_radius, min(720 - ball_radius, player2_pos[1]))
   
    

    # flip() the display to put your work on screen
    pygame.display.flip()

    # limits FPS to 60
    # dt is delta time in seconds since last frame, used for framerate-
    # independent physics.
    dt = clock.tick(60) / 1000

pygame.quit()
sys.exit()
