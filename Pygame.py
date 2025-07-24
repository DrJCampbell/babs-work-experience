import pygame
from sys import exit

pygame.init()
screen = pygame.display.set_mode((800,400))
clock = pygame.time.Clock()

test_surface = pygame.Surface((400,200))
test_surface.fill("Red")

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            exit()

    screen.blit(test_surface,(0,0))

    pygame.display.update()
    clock.tick(60)
    
