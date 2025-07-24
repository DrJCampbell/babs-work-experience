Pygame helps to display images, animations and sounds
through creating a new window, and customizing it
It can also help checking for player input because
(input()) is useless.

(If "p" is before the "." it means "pygame", e.g: pygame.init() )

from ... import ... - Used for importing certain methods (print()) from a module

p.init() - Necessary for starting any pygame code

p.display.set_mode((width,height)) - For creating the window, inside of
                                    the () is arguments you put in, (Put in a
                                     variable)

p.Surface((width,height)) - For putting something on the screen {1}

p.display.update() - Put at the end inside of a while loop, updating the
                    information in the window

screen.blit(surface,(position)) - (blit = block image transfer), for putting 
                                 surfaces on the window screen, in "position" are
                                 coordinates e.g (0,0)

"surfacename".fill("Colour") - "surfacename" is the variable you made {2}, for
                               "Colour", you can put Red or Blue

p.event.get() - Used for getting all the events that the player may do,
               events as in actions (clicking, dragging cursor), use with a for
               loop to check if the user has done any of the actions

p.quit() - Closes the window

p.display.set_caption("name") - Used for setting a game title, in the name space

p.time.Clock() - Used for creating a clock object, used with something else to
                 control the frame rate (step 1), (Stored in a variable)

clock.tick("number") - Number tells the loop to not run faster than that per
                       every second (step 2)

.type - Checking the type of the variable (type of event)

List of event - [QUIT, ACTIVEEVENT, KEYDOWN, KEYUP, MOUSEMOTION, MOUSEBUTTONUP,
                 MOUSEBUTTONDOWN, JOYAXISMOTION, JOYBALLMOTION, JOYHATMOTION,
                 JOYBUTTONUP, JOYBUTTONDOWN, VIDEORESIZE, VIDEOEXPOSE, USEREVENT]

exit() - Imported from the "sys" module, closes all code

p.image.load("Path of your image") - For importing an image, keep the image and
                                     the python file in the same folder




 
                                     
