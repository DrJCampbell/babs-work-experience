import pygame
import sys
import math
import random
from collections import deque

# Initialize Pygame
pygame.init()

# Screen setup
WIDTH, HEIGHT = 1280, 720
CELL_SIZE = 40
cols = WIDTH // CELL_SIZE
rows = HEIGHT // CELL_SIZE
if cols % 2 == 0: cols -= 1
if rows % 2 == 0: rows -= 1

screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Maze Survival")

# Colors
WHITE = (255, 255, 255)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
BLACK = (0, 0, 0)
PURPLE = (128, 0, 128)

# Ball properties
ball_radius = 10
speed = 5
ai_step = 4  # AI movement speed

# Game duration
game_duration = 30000  # 30 seconds
ai_delay = 3000        # Delay before AI starts chasing (in milliseconds)

# Maze generation
def add_extra_openings(maze, chance=0.1):
    for row in range(1, len(maze) - 1):
        for col in range(1, len(maze[0]) - 1):
            if maze[row][col] == 1 and random.random() < chance:
                if (maze[row-1][col] == 0 and maze[row+1][col] == 0) or \
                   (maze[row][col-1] == 0 and maze[row][col+1] == 0):
                    maze[row][col] = 0

def generate_maze(rows, cols):
    maze = [[1 for _ in range(cols)] for _ in range(rows)]

    def carve(x, y):
        directions = [(0, -2), (0, 2), (-2, 0), (2, 0)]
        random.shuffle(directions)
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 1 <= nx < rows - 1 and 1 <= ny < cols - 1 and maze[nx][ny] == 1:
                maze[x + dx // 2][y + dy // 2] = 0
                maze[nx][ny] = 0
                carve(nx, ny)

    maze[1][1] = 0
    carve(1, 1)
    add_extra_openings(maze, chance=0.1)
    return maze

maze = generate_maze(rows, cols)

# Starting positions
ball1_pos = [CELL_SIZE * 1 + CELL_SIZE // 2, CELL_SIZE * 1 + CELL_SIZE // 2]  # Top-left
ball2_pos = [CELL_SIZE * 1 + CELL_SIZE // 2, CELL_SIZE * (rows - 2) + CELL_SIZE // 2]  # Bottom-left
ai_pos = [CELL_SIZE * (cols - 2) + CELL_SIZE // 2, CELL_SIZE * (rows - 2) + CELL_SIZE // 2]  # Bottom-right
ai_target = ai_pos[:]

# Collision detection
def is_wall(pos):
    offsets = [(-ball_radius, -ball_radius), (ball_radius, -ball_radius),
               (-ball_radius, ball_radius), (ball_radius, ball_radius)]
    for dx, dy in offsets:
        x = pos[0] + dx
        y = pos[1] + dy
        col = int(x // CELL_SIZE)
        row = int(y // CELL_SIZE)
        if 0 <= row < len(maze) and 0 <= col < len(maze[0]):
            if maze[row][col] == 1:
                return True
        else:
            return True
    return False

def balls_collide(pos1, pos2, radius):
    return math.hypot(pos1[0] - pos2[0], pos1[1] - pos2[1]) < radius * 2

# BFS pathfinding
def bfs(start, goal, maze):
    start_cell = (int(start[1] // CELL_SIZE), int(start[0] // CELL_SIZE))
    goal_cell = (int(goal[1] // CELL_SIZE), int(goal[0] // CELL_SIZE))
    queue = deque([start_cell])
    came_from = {start_cell: None}

    while queue:
        current = queue.popleft()
        if current == goal_cell:
            break
        for dx, dy in [(-1,0),(1,0),(0,-1),(0,1)]:
            nx = int(current[0] + dy)
            ny = int(current[1] + dx)
            if 0 <= nx < rows and 0 <= ny < cols and maze[nx][ny] == 0:
                neighbor = (nx, ny)
                if neighbor not in came_from:
                    queue.append(neighbor)
                    came_from[neighbor] = current

    path = []
    current = goal_cell
    while current and current in came_from:
        path.append(current)
        current = came_from[current]
    path.reverse()
    return path

def get_next_cell(ai_pos, target_pos, maze):
    path = bfs(ai_pos, target_pos, maze)
    if len(path) > 1:
        next_cell = path[1]
        return [next_cell[1] * CELL_SIZE + CELL_SIZE // 2, next_cell[0] * CELL_SIZE + CELL_SIZE // 2]
    return ai_pos

def smooth_move(current, target, step=2):
    dx = target[0] - current[0]
    dy = target[1] - current[1]
    dist = math.hypot(dx, dy)
    if dist < step:
        return target
    return [current[0] + step * dx / dist, current[1] + step * dy / dist]

# Game loop
clock = pygame.time.Clock()
font = pygame.font.SysFont(None, 48)
running = True
timer_started = False
start_time = 0

while running:
    clock.tick(60)
    screen.fill(WHITE)

    # Draw maze
    for row_idx, row in enumerate(maze):
        for col_idx, cell in enumerate(row):
            if cell == 1:
                pygame.draw.rect(screen, BLACK, (col_idx * CELL_SIZE, row_idx * CELL_SIZE, CELL_SIZE, CELL_SIZE))

    # Start timer after first frame
    if not timer_started:
        start_time = pygame.time.get_ticks()
        timer_started = True

    # Timer
    elapsed = pygame.time.get_ticks() - start_time
    remaining = max(0, game_duration - elapsed)
    timer_text = font.render(f"Time Left: {remaining // 1000}", True, RED)
    screen.blit(timer_text, (20, 20))

    # Event handling
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    keys = pygame.key.get_pressed()

    # Ball 1 movement (WASD)
    new_pos1 = ball1_pos[:]
    if keys[pygame.K_w]: new_pos1[1] -= speed
    if keys[pygame.K_s]: new_pos1[1] += speed
    if keys[pygame.K_a]: new_pos1[0] -= speed
    if keys[pygame.K_d]: new_pos1[0] += speed
    if not is_wall(new_pos1): ball1_pos = new_pos1

    # Ball 2 movement (Arrow keys)
    new_pos2 = ball2_pos[:]
    if keys[pygame.K_UP]: new_pos2[1] -= speed
    if keys[pygame.K_DOWN]: new_pos2[1] += speed
    if keys[pygame.K_LEFT]: new_pos2[0] -= speed
    if keys[pygame.K_RIGHT]: new_pos2[0] += speed
    if not is_wall(new_pos2): ball2_pos = new_pos2

    # AI movement (after delay)
    if elapsed > ai_delay:
        dist1 = math.hypot(ai_pos[0] - ball1_pos[0], ai_pos[1] - ball1_pos[1])
        dist2 = math.hypot(ai_pos[0] - ball2_pos[0], ai_pos[1] - ball2_pos[1])
        target_pos = ball1_pos if dist1 < dist2 else ball2_pos
        ai_target = get_next_cell(ai_pos, target_pos, maze)
        ai_pos = smooth_move(ai_pos, ai_target, step=ai_step)

    # Game over check
    if balls_collide(ball1_pos, ai_pos, ball_radius) or balls_collide(ball2_pos, ai_pos, ball_radius):
        print("Game Over!")
        game_over_text = font.render("Game Over!", True, RED)
        screen.blit(game_over_text, (WIDTH // 2 - 100, HEIGHT // 2))
        pygame.display.flip()
        pygame.time.wait(2000)
        running = False

    # Win condition
    if remaining == 0:
        print("You Win!")
        win_text = font.render("You Win!", True, BLUE)
        screen.blit(win_text, (WIDTH // 2 - 100, HEIGHT // 2))
        pygame.display.flip()
        pygame.time.wait(2000)
        running = False

    # Draw balls (MUST be inside the loop)
    pygame.draw.circle(screen, RED, ball1_pos, ball_radius)
    pygame.draw.circle(screen, BLUE, ball2_pos, ball_radius)
    pygame.draw.circle(screen, PURPLE, ai_pos, ball_radius)

    pygame.display.flip()

pygame.quit()
sys.exit()
