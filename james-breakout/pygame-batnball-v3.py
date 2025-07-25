# Breakout
import pygame
import pygame.font
import random
import time

'''
Ideas with Alex
- unbreakable bricks to create obstacles (or that take multiple)
  hits to break
- exploding super bricks that either take out surrounding bricks
  or send out more balls
- Blackhole bricks - transport ball to a random location?
'''

# pygame setup
pygame.init()
screen = pygame.display.set_mode((1440, 720))
clock = pygame.time.Clock()
dt = clock.tick(60) / 1000
bat_mv_amt = 600
ball_mv_amt_x = 200 # values between 200 and 400 make sense
ball_mv_amt_y = 200 # values between 200 and 400 make sense
ball_mv_incr = 3 # an increase from 200 to 400 in 60 bounces is 3.3
playing = True
score = 0
sample_frame_rate = 500
frame_count = 0
player_last_pos = pygame.Vector2(screen.get_width() / 2, screen.get_height() * 0.9)
player_pos = pygame.Vector2(screen.get_width() / 2, screen.get_height() * 0.9)
ball_pos = pygame.Vector2((screen.get_width()-random.randint(-10,10)) / 2, (screen.get_height()-random.randint(-10,10)) / 2)
ball_vel = pygame.Vector2(ball_mv_amt_x * dt, -ball_mv_amt_y * dt)

def create_bricks(screen, brick_height=30):
    brick_width = int(screen.get_width() / 10)
    #brick_height = 30
    bricks = dict()
    brick_colours = dict()
    brick_counter = 0
    for i in range(0, 90, brick_height):
        for j in range(0, screen.get_width(), brick_width):
            bricks[f"{i}_{j}"] = pygame.Rect(j, i, brick_width, brick_height)
            brick_counter += 1
            #brick_colours[f"{i}_{j}"] = list(["#F5E753", "#49958B", "#40135A"])[brick_counter % 3]
            brick_colours[f"{i}_{j}"] = list(["#5383EC", "#58A55D", "#F2BF42", "#58A55D", "#D85040", "#58A55D", "#F2BF42"])[brick_counter % 7]
    return(bricks, brick_colours)

bricks, brick_colours = create_bricks(screen)


# set up on-screen text
font = pygame.font.Font(None, 64)
text_oob = font.render("Out of bounds", True, (10, 10, 10))
textpos_oob = text_oob.get_rect(centerx=screen.get_width() / 2, y=screen.get_height() / 2)
font_l_text = font.render("left", True, (255, 255, 255))
font_r_text = font.render("right", True, (255, 255, 255))
font_lr_pos = font_l_text.get_rect(centerx=screen.get_width() / 2, y=screen.get_height() / 2)
while playing:
    # poll for events
    # pygame.QUIT event means the user clicked X to close your window
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            playing = False
    
    frame_count += 1

    # fill the screen with a color to wipe away anything from last frame
    screen.fill("black")

    for brick in bricks:
        pygame.draw.rect(screen, brick_colours[brick], bricks[brick])


    bat_rect = pygame.Rect(player_pos.x, player_pos.y, 100, 10)
    pygame.draw.rect(screen, "white", bat_rect)
    
    pygame.draw.circle(screen, "white", ball_pos, 10)

    if bat_rect.collidepoint(ball_pos):
        ball_mv_amt_x += ball_mv_incr
        ball_mv_amt_y += ball_mv_incr
        
        # check if bat is going left, right or stationary
        # if going in opposite direction to ball:
        # - decrease ball_mv_amt_y
        # - increase ball_mv_amt_x
        # reverse these for bat in the same direction

        #screen.blit(font_l_text, font_lr_pos)
        # bat moving left
        if(player_last_pos.x > player_pos.x):
            ball_mv_amt_y *= 0.5
            ball_mv_amt_x *= 2
            screen.blit(font_l_text, font_lr_pos)
            #time.sleep(1)
        elif(player_last_pos.x < player_pos.x):
            ball_mv_amt_y *= 2
            ball_mv_amt_x *= 0.5
            screen.blit(font_r_text, font_lr_pos)
            #time.sleep(1)
        # else:
        #     ball_mv_amt_y = 0 # 2
        #     ball_mv_amt_x = 0 # 0.5
        #     font_pos_info = font.render(f"prev: {player_last_pos.x} curr: {player_pos.x}", True, (255, 255, 255))
        #     font_pos_pos = font_pos_info.get_rect(centerx=screen.get_width() / 2, y=screen.get_height() / 2)
        #     screen.blit(font_pos_info, font_pos_pos)

        ball_vel.y = -ball_mv_amt_y * dt
    
    if ball_pos.x <= 0:
        ball_vel.x = ball_mv_amt_x * dt
    if ball_pos.x >= screen.get_width():
        ball_vel.x = -ball_mv_amt_x * dt

    if ball_pos.y >= screen.get_height():

        # Game over - show menu
        waiting = True
        while waiting:
            # show menu stuff
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    waiting = False
                    playing = False
            
            font = pygame.font.Font(None, 64)
            text_menu = font.render("Press P to play or Q to quit", True, (255, 255, 255))
            textpos_menu = text_menu.get_rect(centerx=screen.get_width() / 2, y=screen.get_height() / 2)
            screen.blit(text_menu, textpos_menu)

            keys = pygame.key.get_pressed()
            if keys[pygame.K_q]:
                waiting = False
                playing = False
            if keys[pygame.K_p]:
                score = 0
                bricks, brick_colours = create_bricks(screen)
                ball_pos = pygame.Vector2((screen.get_width()-random.randint(-10,10)) / 2, (screen.get_height()-random.randint(-10,10)) / 2)
                ball_mv_amt_x = 200
                ball_mv_amt_y = 200
                ball_vel = pygame.Vector2(ball_mv_amt_x * dt, -ball_mv_amt_y * dt)
                waiting = False
            
            pygame.display.flip()

    if ball_pos.y <= 0:
        ball_vel.y = ball_mv_amt_y * dt
    brick_to_remove = None
    for brick in bricks:
        if bricks[brick].collidepoint(ball_pos):
            ball_vel.y = -ball_vel.y
            brick_to_remove = brick
    if brick_to_remove:
        del(bricks[brick_to_remove])
        score += 1
    
    # we determine the trajectory of the ball here
    ball_pos.x = ball_pos.x + ball_vel.x
    ball_pos.y = ball_pos.y + ball_vel.y

    if(frame_count % sample_frame_rate) == 0:
        player_last_pos = player_pos
    
    #font_pos_info = font.render(f"prev: {player_last_pos.x}\\ncurr: {player_pos.x}", True, (255, 255, 255))
    #font_pos_pos = font_pos_info.get_rect(centerx=screen.get_width() / 2, y=screen.get_height() / 2)
    #screen.blit(font_pos_info, font_pos_pos)


    keys = pygame.key.get_pressed()
    if keys[pygame.K_a] or keys[pygame.K_LEFT]:
        if player_pos.x > 0:
            player_pos.x -= bat_mv_amt * dt
    if keys[pygame.K_d] or keys[pygame.K_RIGHT]:
        if (player_pos.x + 100) < screen.get_width():
            player_pos.x += bat_mv_amt * dt
    
    

    if len(bricks) == 0:
        bricks, brick_colours = create_bricks(screen)
    
    text_score = font.render(f"score: {score}", True, (255, 255, 255))
    textpos_score = text_score.get_rect(centerx=screen.get_width() * 0.80, y=screen.get_height() * 0.90)
    screen.blit(text_score, textpos_score)

    # flip() the display to put your work on screen
    pygame.display.flip()

    # limits FPS to 60
    # dt is delta time in seconds since last frame, used for framerate-
    # independent physics.
    dt = clock.tick(60) / 1000


pygame.quit()

