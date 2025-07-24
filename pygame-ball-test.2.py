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
pygame.display.set_caption("Chase the AI")

# Colors
WHITE = (255, 255, 255)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
BLACK = (0, 0, 0)
PURPLE = (128, 0, 128)

# Ball properties
ball_radius = 10
speed = 5
ai_step = 4.5

# Game duration
game_duration = 30000  # 30 seconds
ai_delay = 500         # Delay before AI starts moving

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
ball1_pos = [CELL_SIZE * 1 + CELL_SIZE // 2, CELL_SIZE * 1 + CELL_SIZE // 2]
ball2_pos = [CELL_SIZE * (cols - 2) + CELL_SIZE // 2, CELL_SIZE * 1 + CELL_SIZE // 2]
ai_pos = [CELL_SIZE * (cols // 2) + CELL_SIZE // 2, CELL_SIZE * (rows - 2) + CELL_SIZE // 2]

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

ai_move_timer = 0
ai_move_interval = 200  # AI picks new target every 200ms

# AI path tracking
ai_path = []
ai_path_index = 0

while running:
    dt = clock.tick(60)
    screen.fill(WHITE)

    # Draw maze
    for row_idx, row in enumerate(maze):
        for col_idx, cell in enumerate(row):
            if cell == 1:
                pygame.draw.rect(screen, BLACK, (col_idx * CELL_SIZE, row_idx * CELL_SIZE, CELL_SIZE, CELL_SIZE))

    if not timer_started:
        start_time = pygame.time.get_ticks()
        timer_started = True

    elapsed = pygame.time.get_ticks() - start_time
    remaining = max(0, game_duration - elapsed)
    timer_text = font.render(f"Time Left: {remaining // 1000}", True, RED)
    screen.blit(timer_text, (20, 20))

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    keys = pygame.key.get_pressed()

    # Player 1 movement (WASD)
    new_pos1 = ball1_pos[:]
    if keys[pygame.K_w]: new_pos1[1] -= speed
    if keys[pygame.K_s]: new_pos1[1] += speed
    if keys[pygame.K_a]: new_pos1[0] -= speed
    if keys[pygame.K_d]: new_pos1[0] += speed
    if not is_wall(new_pos1): ball1_pos = new_pos1

    # Player 2 movement (Arrow keys)
    new_pos2 = ball2_pos[:]
    if keys[pygame.K_UP]: new_pos2[1] -= speed
    if keys[pygame.K_DOWN]: new_pos2[1] += speed
    if keys[pygame.K_LEFT]: new_pos2[0] -= speed
    if keys[pygame.K_RIGHT]: new_pos2[0] += speed
    if not is_wall(new_pos2): ball2_pos = new_pos2

    # AI movement (after delay)
    if elapsed > ai_delay:
        ai_move_timer += dt
        if ai_move_timer >= ai_move_interval or not ai_path:
            ai_move_timer = 0
            # Choose best target maximizing distance from players
            best_target = (int(ai_pos[0] // CELL_SIZE), int(ai_pos[1] // CELL_SIZE))
            max_dist = -1
            best_path = []

            for _ in range(100):
                row = random.randint(1, rows - 2)
                col = random.randint(1, cols - 2)
                if maze[row][col] == 0:
                    cell_center = [col * CELL_SIZE + CELL_SIZE // 2, row * CELL_SIZE + CELL_SIZE // 2]
                    dist1 = math.hypot(cell_center[0] - ball1_pos[0], cell_center[1] - ball1_pos[1])
                    dist2 = math.hypot(cell_center[0] - ball2_pos[0], cell_center[1] - ball2_pos[1])
                    closest_dist = min(dist1, dist2)
                    path = bfs(ai_pos, cell_center, maze)
                    if path and closest_dist > max_dist:
                        max_dist = closest_dist
                        best_target = (row, col)
                        best_path = path

            ai_path = best_path
            ai_path_index = 0

        # Move AI along the path if it exists
        if ai_path and ai_path_index < len(ai_path):
            # Convert path cell to actual pixel position
            next_cell = ai_path[ai_path_index]
            next_pos = [next_cell[1] * CELL_SIZE + CELL_SIZE // 2, next_cell[0] * CELL_SIZE + CELL_SIZE // 2]
            ai_pos = smooth_move(ai_pos, next_pos, ai_step)
            # If AI is close enough to next cell, go to next path index
            if math.hypot(ai_pos[0] - next_pos[0], ai_pos[1] - next_pos[1]) < 1.5:
                ai_path_index += 1
        else:
            # No path or reached end, stay put until next target chosen
            pass

    # Check catch condition
    if balls_collide(ball1_pos, ai_pos, ball_radius) or balls_collide(ball2_pos, ai_pos, ball_radius):
        print("You Win!")
        win_text = font.render("You Win!", True, BLUE)
        screen.blit(win_text, (WIDTH // 2 - 100, HEIGHT // 2))
        pygame.display.flip()
        pygame.time.wait(2000)
        running = False

    if remaining == 0:
        print("Game Over!")
        game_over_text = font.render("Game Over!", True, RED)
        screen.blit(game_over_text, (WIDTH // 2 - 100, HEIGHT // 2))
        pygame.display.flip()
        pygame.time.wait(2000)
        running = False

    # Draw balls
    pygame.draw.circle(screen, RED, [int(x) for x in ball1_pos], ball_radius)
    pygame.draw.circle(screen, BLUE, [int(x) for x in ball2_pos], ball_radius)
    pygame.draw.circle(screen, PURPLE, [int(x) for x in ai_pos], ball_radius)

    pygame.display.flip()

pygame.quit()
sys.exit()
