from moviepy.editor import VideoFileClip

def mp4_to_gif(mp4_file, gif_file):
    # Load the video file
    video = VideoFileClip(mp4_file)
    
    # Optional: Resize the video to reduce dimensions (you can change the scale factor as needed)
    video = video.resize(height=300, width=300)  # resize height to 320px, width is adjusted proportionally
    
    # Optional: Set a frame rate to downsample (lowering this will reduce the GIF file size)
    video = video.set_fps(30)  # reduce the frame rate (default is 24 or higher)
    
    # Write the gif file
    video.write_gif(gif_file)

# Example usage
mp4_to_gif("Tiny FB MPC.mp4", "Tiny FB MPC.gif")
