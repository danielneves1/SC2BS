

from PIL import Image, ImageSequence
from tqdm import tqdm




def concatenate_gifs(gif_paths, output_path, duration=1, loop=0):
    """
    Concatenate multiple GIFs one after another and save as a new animated GIF.

    Parameters:
    - gif_paths: List of paths to input GIF files.
    - output_path: Output path for the combined GIF.
    - duration: Duration between frames in milliseconds (default is 500ms).
    - loop: Number of loops (0 means infinite).
    """
    frames = []

    # Loop through each GIF and extract its frames
    for gif_path in tqdm(gif_paths):
        gif = Image.open(gif_path)

        # Extract each frame from the GIF and add it to the frame list

        for frame in ImageSequence.Iterator(gif):
            # Ensure that we have a full copy of the frame before appending
            frames.append(frame.copy())

    # Save the concatenated frames as a single GIF
    frames[0].save(output_path,
                   save_all=True,
                   append_images=frames[1:], loop=loop) # All frames after the first one
                #   duration=100,  # Duration between frames
                     # Looping (0 means infinite loop)


if __name__=="__main__":


    #example

    workers = 2
    runs = 2

    n_gifs = int(runs * workers) # we will combine n_gifs
    giflocation = "/Output/Gifs/norm/"
    gifname = "u_" #in this example each gif is labelled as u_0.gif, u_1.gif, u_2.gif etc...
    gif_paths = [giflocation + gifname + str(s)+".gif" for s in range(0, n_gifs)]# Paths to the gifs you wish to join

    output_path =  giflocation + "u.gif"  # Path to save the combined gif
    concatenate_gifs(gif_paths, output_path)

    giflocation = "/Output/Gifs/sqrd/"
    gifname = "u2_"  # in this example each gif is labelled as u2_0.gif, u2_1.gif, u2_2.gif etc...
    gif_paths = [giflocation + gifname + str(s) + ".gif" for s in
                 range(0, n_gifs)]  # Paths to the gifs you wish to join

    output_path =  giflocation  + "u2.gif"  # Path to save the combined gif
    concatenate_gifs(gif_paths, output_path)
